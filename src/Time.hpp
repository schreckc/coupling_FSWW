#ifndef TIMES_HPP
#define TIMES_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
#include <time.h>
using namespace boost::posix_time;

class Times {
  
public :
  static Times *TIMES;
  static Times *TIMES_UP;

  enum timing_t {deform_time_ = 0, 
		 solve_time_,
		 nTimes};

private :
  boost::posix_time::ptime init_time[nTimes];
  double loop_time[nTimes];
  double time_sum[nTimes];

  //count
  unsigned int frame;
  unsigned int step;

public :
  Times() {
  init();
}

  
  void init() {
  for (unsigned int i = 0; i < nTimes; ++i) {
    init_time[i] = not_a_date_time;
    loop_time[i] = 0;
    time_sum[i] = 0;
  }
  frame = 0;
  step = 0;
}

  void tick(unsigned int i) {
  init_time[i] = microsec_clock::local_time();
}

  void tock(unsigned int i) {
  assert(init_time[i] != not_a_date_time);
  ptime t_end = microsec_clock::local_time();
   loop_time[i] += (t_end - init_time[i]).total_microseconds()*1e-6;
   init_time[i] = not_a_date_time;
}

  double getTime(unsigned int i) {
  return loop_time[i];
}

  double getAverageTime(unsigned int i) {
  //  INFO("time_sum[i]", time_sum[i]<<" "<<cpt<<" "<<time_sum[i]/((double)cpt));
  return time_sum[i]/((double)frame);
}

  double getAverageTimeByStep(unsigned int i) {
  //  INFO("time_sum[i]", time_sum[i]<<" "<<cpt<<" "<<time_sum[i]/((double)cpt));
  return time_sum[i]/((double)step);
}

  double getAverageNbStepsByFrame() {
  return (double) step / ((double) frame);
}

  
  void next_loop() {
  for (unsigned int i = 0; i < nTimes; ++i) {
    time_sum[i] += loop_time[i];
    init_time[i] = not_a_date_time;
    loop_time[i] = 0;
  }

  ++frame;
}
  void next_step()  {
  ++step;
}

};

#endif //TIMES_HPP
