#include "async_output.h"

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace {
std::unique_ptr<PARTICLE> makeParticleSnapshot(const PARTICLE& source) {
    auto snapshot = std::make_unique<PARTICLE>();
    const int count = source.Num;
    snapshot->Num = count;
    snapshot->Mass = new double[count];
    snapshot->Radius = new double[count];
    snapshot->Number = new int[count];
    snapshot->Status = new int[count];
    snapshot->Pos = new vector3d[count];
    snapshot->Vel = new vector3d[count];
    snapshot->AngSpd = new vector3d[count];
    snapshot->Quat = new vector4d[count];
    std::copy(source.Mass, source.Mass + count, snapshot->Mass);
    std::copy(source.Radius, source.Radius + count, snapshot->Radius);
    std::copy(source.Number, source.Number + count, snapshot->Number);
    std::copy(source.Status, source.Status + count, snapshot->Status);
    snapshot->maxRadius = source.maxRadius;
    snapshot->MeshSize = source.MeshSize;
    return snapshot;
}
enum class FrameState { Free, Filling, Queued, Writing };
}

OutputFrame::OutputFrame(const PARTICLE& source, int bodyCount,
                         const AsyncOutputOptions& options)
    : options_(options),
      particles_(options.particleState ? makeParticleSnapshot(source) : nullptr),
      forces_((options.particleForce || (options.bodyForce && bodyCount > 0))
                  ? std::make_unique<FORCE>(source.Num, bodyCount) : nullptr),
      bodyStorage_(options.bodyState && bodyCount > 0
                       ? std::make_unique<BODY[]>(bodyCount) : nullptr),
      bodyForces_(options.bodyState && bodyCount > 0 ? bodyCount : 0),
      bodyForceRows_(options.bodyState && bodyCount > 0 ? bodyCount : 0),
      bodyImpulses_(options.bodyState && bodyCount > 0 ? bodyCount : 0),
      bodyImpulseRows_(options.bodyState && bodyCount > 0 ? bodyCount : 0),
      bodyCount_(bodyCount) {
    bodies_.Num = options.bodyState ? bodyCount : 0;
    bodies_.body = bodyStorage_.get();
    for (std::size_t i = 0; i < bodyForceRows_.size(); ++i) {
        bodyForceRows_[i] = bodyForces_[i].data();
        bodyImpulseRows_[i] = bodyImpulses_[i].data();
    }
}
OutputFrame::~OutputFrame() { bodies_.body = nullptr; bodies_.Num = 0; }
PARTICLE& OutputFrame::particles() {
    if (!particles_) throw std::logic_error("particle state output is disabled");
    return *particles_;
}
FORCE& OutputFrame::forces() {
    if (!forces_) throw std::logic_error("force output is disabled");
    return *forces_;
}
void OutputFrame::captureBodies(
    const BODYSET& source, double** loads, double** impulses
) {
    if (source.Num != bodyCount_)
        throw std::runtime_error("body count changed while capturing output");
    if (!options_.bodyState) return;
    for (int i = 0; i < bodyCount_; ++i) {
        bodyStorage_[i].MassCenter = source.body[i].MassCenter;
        bodyStorage_[i].Vel = source.body[i].Vel;
        bodyStorage_[i].AngularVel = source.body[i].AngularVel;
        bodyStorage_[i].orien = source.body[i].orien;
        for (int k = 0; k < 6; ++k) {
            bodyForces_[i][k] = loads[i][k];
            bodyImpulses_[i][k] = impulses[i][k];
        }
    }
}
void OutputFrame::setMetadata(int step, double time) { step_ = step; time_ = time; }
void OutputFrame::writeFiles() {
    if (options_.particleForce) {
        if (bodyCount_ > 0) forces_->ParticlesForceOutput(step_, forces_->ptNum, bodyCount_);
        else forces_->ParticlesForceOutput(step_, forces_->ptNum);
    }
    if (options_.bodyForce && bodyCount_ > 0)
        forces_->BodysetForceOutput(step_, forces_->ptNum, bodyCount_);
    if (options_.particleState) particles_->StateOutput(step_);
    if (options_.bodyState && bodyCount_ > 0)
        bodies_.StateOutput(
            time_, bodyForceRows_.data(), step_, bodyImpulseRows_.data());
}

class AsyncOutputWriter::Impl {
public:
    Impl(const PARTICLE& particles, int bodyCount, const AsyncOutputOptions& options) {
        for (std::size_t i=0;i<frames_.size();++i) {
            frames_[i]=std::make_unique<OutputFrame>(particles,bodyCount,options);
            states_[i]=FrameState::Free;
        }
        worker_=std::thread(&Impl::workerLoop,this);
    }
    ~Impl() {
        { std::lock_guard<std::mutex> lock(mutex_); stop_=true; }
        queued_.notify_all();
        if(worker_.joinable()) worker_.join();
    }
    OutputFrame& acquire() {
        std::unique_lock<std::mutex> lock(mutex_);
        free_.wait(lock,[this]{return error_||findFree()<frames_.size();});
        rethrowError();
        const auto index=findFree(); states_[index]=FrameState::Filling;
        return *frames_[index];
    }
    void submit(OutputFrame& frame) {
        std::lock_guard<std::mutex> lock(mutex_); rethrowError();
        const auto index=findFrame(frame);
        if(index==frames_.size()||states_[index]!=FrameState::Filling)
            throw std::logic_error("invalid asynchronous output frame submission");
        states_[index]=FrameState::Queued; queue_.push_back(index); queued_.notify_one();
    }
    void finish() {
        std::unique_lock<std::mutex> lock(mutex_);
        idle_.wait(lock,[this]{return error_||(queue_.empty()&&writes_==0);});
        rethrowError();
    }
private:
    std::size_t findFree() const {
        for(std::size_t i=0;i<states_.size();++i) if(states_[i]==FrameState::Free) return i;
        return states_.size();
    }
    std::size_t findFrame(const OutputFrame& frame) const {
        for(std::size_t i=0;i<frames_.size();++i) if(frames_[i].get()==&frame) return i;
        return frames_.size();
    }
    void rethrowError() const { if(error_) std::rethrow_exception(error_); }
    void workerLoop() {
        while(true) {
            std::size_t index=0;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                queued_.wait(lock,[this]{return stop_||error_||!queue_.empty();});
                if(error_) return;
                if(queue_.empty()){if(stop_)return;continue;}
                index=queue_.front();queue_.pop_front();states_[index]=FrameState::Writing;++writes_;
            }
            try { frames_[index]->writeFiles(); }
            catch(...) {
                std::lock_guard<std::mutex> lock(mutex_);
                error_=std::current_exception();states_[index]=FrameState::Free;--writes_;
                while(!queue_.empty()){states_[queue_.front()]=FrameState::Free;queue_.pop_front();}
                free_.notify_all();idle_.notify_all();return;
            }
            {
                std::lock_guard<std::mutex> lock(mutex_);
                states_[index]=FrameState::Free;--writes_;free_.notify_one();
                if(queue_.empty()&&writes_==0)idle_.notify_all();
            }
        }
    }
    std::array<std::unique_ptr<OutputFrame>,2> frames_;
    std::array<FrameState,2> states_{};
    std::deque<std::size_t> queue_;
    std::thread worker_;
    std::mutex mutex_;
    std::condition_variable free_,queued_,idle_;
    std::exception_ptr error_;
    std::size_t writes_=0;
    bool stop_=false;
};

AsyncOutputWriter::AsyncOutputWriter(const PARTICLE& p,int n,const AsyncOutputOptions& o)
    :impl_(std::make_unique<Impl>(p,n,o)){}
AsyncOutputWriter::~AsyncOutputWriter()=default;
OutputFrame& AsyncOutputWriter::acquire(){return impl_->acquire();}
void AsyncOutputWriter::submit(OutputFrame& f){impl_->submit(f);}
void AsyncOutputWriter::finish(){impl_->finish();}
