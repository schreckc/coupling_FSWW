#include "definitions.hpp"
#include "InputPoint.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <vector>
#include <Eigen/SVD>

int main(int argc, char **argv) {
  std::string file_name;
  std::string merge_name;
  int start = 0;
  int stop = 2000;
  bool plot_ = false;
  bool acc_max = false;
  bool merge_ = false;
  int span = 6;
  bool follow_ = false;
  bool fft_ = false;

  float endx = 0, beginx = 0;
  float ymin = 0, ymax = 0;
  
  for (int i = 1;  i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "-file") {
      assert(argc > i+1);
      file_name = argv[i+1];
      ++i;
    }
    if (s == "-start") {
      assert(argc > i+1);
      start = atoi(argv[i+1]);
      ++i;
    }
    if (s == "-stop") {
      assert(argc > i+1);
      stop = atoi(argv[i+1]);
      ++i;
    }
    if (s == "-plot") {
      std::cout<<"PLOT"<<std::endl;
      plot_ = true;
    }
    if (s == "-acc_max") {
      acc_max = true;
    }
    if (s == "-merge") {
      assert(argc > i+1);
      merge_name = argv[i+1];
      merge_ = true;
      ++i;
    }
    if (s == "-span") {
      assert(argc > i+1);
      span = atoi(argv[i+1]);
      ++i;
    }
    if (s == "-follow") {
      follow_ = true;
    }
    if (s == "-fft") {
      fft_ = true;
    }
  }

  if (follow_) {
    std::vector<std::list<float> > max_h(stop-start);
    std::vector<std::list<float> > max_x(stop-start);

    std::list<std::list<float> > maxf_h;
    std::list<std::list<float> > maxf_x;
    int m = 0;
  
    for (uint ind = start; ind < stop; ++ind) {
      std::ifstream  file;
      std::list<float> x_list;
      std::list<float> h_list;
    
      std::stringstream ss; 
      ss <<file_name<<ind<<".txt";
      std::string str(ss.str());
      file.open(str);
    
      if (!file.good()) {
	std::cout<<"FILE NOT GOOD: -- "<<str<<std::endl;
	file.close();
	//      return -1;
      } else {
	int nb_pts = 0;
	int nb_max = 0;
	int nb_min = 0;
	while (!file.eof()) {
	  float x, h;
	  file >> x >> h;
	  x_list.push_front(x);
	  h_list.push_front(h);
	  if (h > ymax) {
	    ymax = h;
	    m = ind;
	  }
	  if (h < ymin) {
	    ymin = h;
	  }
	  // std::cout<<x<<" "<<h<<std::endl;
	  ++nb_pts;
	}

	std::vector<float> x_v(nb_pts);
	std::vector<float> h_v(nb_pts);
	std::list<float>::iterator itx = x_list.begin(), ith = h_list.begin();
	for (int j = 0; itx !=  x_list.end(); ++itx, ++ith, ++j) {
	  x_v[nb_pts - 1 - j] = *itx;
	  h_v[nb_pts - 1 - j] = *ith;
	  //       std::cout<<j<<" "<<x_v[j]<<" "<<h_v[j]<<std::endl;
	}
	if (endx == 0) {
	  endx = x_v[0];
	  beginx = x_v[nb_pts-1];
	}
	for (int j = span; j < nb_pts - span; ++j) {
	  float cur = h_v[j];
	  int k;
	  //      std::cout<<j<<" "<<cur<<" "<<span<<std::endl;
	  for (k = 1; k <= span; ++k) {
	    if (cur <= h_v[j-k] || cur < h_v[j+k]) {
	      //	std::cout<<"BREAK"<<std::endl;
	      break;
	    }
	  }
	  if (k > span) {
	    max_x[ind-start].push_back(x_v[j]);
	    max_h[ind-start].push_back(cur);
	    ++nb_max;
	  }
	}
	// std::cout<<"NB PTS "<<nb_pts<<" "<<ind<<std::endl;
	// std::cout<<"NB MAXS "<<nb_max<<" "<<max_x[ind-start].size()<<std::endl;
	file.close();
      }
    }
    std::cout<<"FILE WITH MAX "<<m<<" "<<ymax<<std::endl;
    for (int i = start; i < stop; ++i) {
      if (!max_x[i-start].empty()) {
	max_x[i-start].pop_back();
	max_h[i-start].pop_back();
      }
    }
    
    uint count = 0;
    //    std::cout<<"nmax BEGN  "<<max_x[0].size()<<" "<<start<<" "<<stop<<std::endl;      
    for (int i = start; i < stop; ++i) {
      //       std::cout<<"nmax BEGN  "<<i<<" "<<max_x[i-start].size()<<std::endl;      
      // if (max_x[i-start].empty()) {
      // 	break;
      // }
      while (!max_x[i-start].empty()) {
	float prev = max_x[i-start].front();
	std::list<float> nmax_x;
	std::list<float> nmax_h;
	nmax_x.push_back(prev);
	nmax_h.push_back(max_h[i-start].front());
	for (int j = i+1; j < stop; ++j) {
	  // std::list<float> &curx = max_x[j-start];
	  // std::list<float> &curh = max_h[j-start];
	  std::list<float>::iterator itx = max_x[j-start].begin(), ith = max_h[j-start].begin();
	  for (int k = 0; itx !=  max_x[j-start].end(); ++itx, ++ith, ++k) {
	    float d = *itx - prev;
	    if (d >= 0) {
	      break;
	    }
	  }
	  if (itx == max_x[j-start].end()) {
	    break;
	  }
	  nmax_x.push_back(*itx);
	  nmax_h.push_back(*ith);
	  prev = *itx;
	  max_x[j-start].erase(itx);
	  max_h[j-start].erase(ith);
	  //  std::cout<<"nmax "<<nmax_x.size()<<" "<<j<<" "<<max_x[j-start].size()<<std::endl;      
	}
	maxf_x.push_back(nmax_x);
	maxf_h.push_back(nmax_h);
	max_x[i-start].pop_front();
	max_h[i-start].pop_front();
	++count;
	//	std::cout<<"nmax END "<<nmax_x.size()<<" "<<count<<" "<<maxf_x.size()<<" "<<max_x[i-start].size()<<std::endl;      

      }
    }

    std::list<std::list<float> >::iterator itx_l = maxf_x.begin(), ith_l = maxf_h.begin();
    uint i = 0;
    for (; itx_l != maxf_x.end(); ++itx_l, ++ith_l) {
      if ((*itx_l).size() > 200) {
	std::stringstream ss_m;
	ss_m <<file_name<<i<<"_follow.txt";
	std::string str_m(ss_m.str());
	std::ofstream  file_m;
	file_m.open(str_m);
	
	std::list<float> curx = *itx_l;
	std::list<float> curh = *ith_l;
	std::list<float>::iterator itx = curx.begin(), ith = curh.begin();
	// ++itx; ++ith;
	// ++itx; ++ith;
	for (int k = 0; itx !=  curx.end(); ++itx, ++ith, ++k) {
	  file_m << *itx<<" "<<*ith<<"\n";
	}
	file_m.close();
	std::cout<<"nmax "<<curx.size()<<" "<<i<<std::endl;  
	++i;
      }
    }
    ////////////////////////////  
  } else {
  
    std::stringstream ss, ss_max, ss_min; 
    ss_max <<file_name<<"max.txt";
    ss_min <<file_name<<"min.txt";
    std::string str_max(ss_max.str());
    std::string str_min(ss_min.str());
    std::ofstream  file_max;
    std::ofstream  file_min;
   
    file_max.open(str_max);
    file_min.open(str_min);

    //  float ymin = 0, ymax = 0;

    std::list<float> max_x;
    std::list<float> max_h;
    std::list<float> min_x;
    std::list<float> min_h;
    //  float endx = 0, beginx = 0;
    float mean = 0;
    int n_wl = 0;
    InputPoint ip(512, 0.03);
    int m = 0, mm = 0;
    
    float suma = 0;
    float sumb = 0;
    int n_sigma = 0;

    for (uint ind = start; ind < stop; ++ind) {
      if (!acc_max) {
	max_h.clear();
	min_h.clear();
	max_x.clear();
	min_x.clear();
      }
      std::ifstream  file;
      std::list<float> x_list;
      std::list<float> h_list;
  
      std::stringstream ss; 
      ss <<file_name<<ind<<".txt";
      std::string str(ss.str());
      file.open(str);

      //  std::cout<<"read "<<str<<std::endl;
  
      std::ofstream file_t;
      //  if (acc_max) {
      std::stringstream ss_t; 
      ss_t <<file_name<<"max_"<<ind<<".txt";
      std::string str_t(ss_t.str());
      file_t.open(str_t);
      // }

      if (!file.good() || !file_max.good()) {
	std::cout<<"FILE NOT GOOD "<<str<<ind<<std::endl;
	file_t.close();
	file.close();
	//    return -1;
      } else {
	int nb_pts = 0;
	int nb_max = 0;
	int nb_min = 0;

	while (!file.eof()) {
	  float x, h;
	  file >> x >> h;
	  if (file.fail()) {
	    //   std::cout<<"BREAK (file fail) "<<str<<std::endl;
	    break;
	  }
	  x_list.push_back(x);
	  h_list.push_back(h);
	  if (h > ymax) {
	    ymax = h;
	    m = ind;
	  }
	  if (h < ymin) {
	    ymin = h;
	    mm = ind;
	  }
	  // std::cout<<x<<" "<<h<<std::endl;
	  ++nb_pts;
	}
	if (nb_pts == 0) {
	  std::cout<<"BREAK (nb_pts = 0) "<<str<<std::endl;
	  break;
	}
	//  std::cout<<"nb pts "<<nb_pts<<std::endl;
	std::vector<float> x_v(nb_pts);
	std::vector<float> h_v(nb_pts);
	std::list<float>::iterator itx = x_list.begin(), ith = h_list.begin();
	for (int i = 0; /*itx !=  x_list.end()*/i < nb_pts; ++itx, ++ith, ++i) {
	  x_v[i] = *itx;
	  h_v[i] = *ith;
	  // std::cout<<i<<" "<<x_v[i]<<" "<<h_v[i]<<std::endl;
	}
	if (fft_) {
	  if (ind == start) {
	    std::stringstream ss_ip; 
	    ss_ip <<file_name<<"ip";
	    std::string str_ip(ss_ip.str());
	    float xm = x_v[nb_pts/2];
	    ip.setPos(xm, 0);
	    ip.setName(str_ip);
	    ip.rec = true;
	  }
	  ip.update(h_v[nb_pts/2]);
	}
	if (endx == 0) {
	  endx = x_v[0];
	  beginx = x_v[nb_pts-1];
	}
	for (int i = span; i < nb_pts/2 - span; ++i) {
	  float cur = h_v[i];
	  int k;
	  for (k = 1; k <= span; ++k) {
	    if (h_v[i] <= 0 || cur <= h_v[i-k] || cur < h_v[i+k]) {
	      break;
	    }
	  }
	  if (k > span) {
	    max_x.push_back(x_v[i]);
	    max_h.push_back(cur);
	    file_max << x_v[i] <<" "<<cur<<"\n";
     
	    ++nb_max;
	  }
	  for (k = 1; k <= span; ++k) {
	    if (cur >= h_v[i-k] || cur > h_v[i+k]) {
	      break;
	    }
	  }
	  if (k > span) {
	    min_x.push_back(x_v[i]);
	    min_h.push_back(cur);
	    file_min << x_v[i] <<" "<<cur<<"\n";
	    ++nb_min;
	  }
	}
	//  if (acc_max) {
	ith = max_h.begin(); itx = max_x.begin();
	// ++ith; ++itx;
	// ++ith; ++itx;
	for (int k = 0; ith != max_h.end(); ++ith, ++itx, ++k) {
	  file_t << *itx <<" "<<*ith<<" "<<log(*ith)<<"\n";
	  // if (k > nb_max - 2) {
	  // 	break;
	  // }
	  //   }
	}
	file_t.close();
	file.close();
	// std::cout<<"NB PTS "<<nb_pts<<" "<<str<<std::endl;
	// std::cout<<"NB MAXS "<<nb_max<<std::endl;

	ith = max_h.begin(); itx = max_x.begin();
	float prevh = *ith;
	float prevx = *itx;
	float sum = 0;
	int n = 0;
	++ith; ++itx;
	for (; ith != max_h.end(); ++ith, ++itx) {
	  float x = *itx;
	  sum += prevx - x;
	  prevx = x;
	  prevh = *ith;
	  ++n;
	}
	if (n != 0) {
	  float mean_wl = sum / (float)n;
	  mean += mean_wl;
	  std::cout<<"MEAN WL "<<mean_wl<<std::endl;
	  ++n_wl;
	}
      }
      if (acc_max) {
	std::list<float>::iterator itx = max_x.begin(), ith = max_h.begin();
	MatrixXf M(max_h.size(), 2);
	VectorXf b(max_h.size());
	for (int k = 0; ith != max_h.end(); ++ith, ++itx, ++k) {
	  if (*ith > 0) {
	    M(k, 0) = -(*itx);
	    M(k, 1) = 1;
	    b(k) = log(*ith);
	  } else {
	    M(k, 0) = 0;
	    M(k, 1) = 0;
	    b(k) = 0;
	  }
	}
	// std::cout<<M<<std::endl;
	// std::cout<<b<<std::endl;
	Eigen::BDCSVD<MatrixXf> svd(M, ComputeThinU | ComputeThinV);
	VEC2 x;
	x = svd.solve(b);
	std::cout<<"B "<<x(0)<<"\nA "<<exp(x(1))<<std::endl;
	suma += exp(x(1));
	sumb += x(0);
	++n_sigma;
      }
    }
  
    // std::cout<<"FILE WITH MAX "<<m<<" "<<ymax<<std::endl;
    // std::cout<<"FILE WITH MIN "<<mm<<" "<<ymin<<std::endl;
    file_max.close();
    file_min.close();
    mean /= (float)(n_wl);
    std::cout<<"MEAN ALL  "<<mean<<"  A"<<suma/n_sigma<<" B"<<sumb/n_sigma<<std::endl;
  
    if (plot_) {
      std::ofstream stream_plot;
      std::stringstream ss_plot; 
      ss_plot <<file_name<<".plot";
      std::string str_plot(ss_plot.str());
      std::cout<<"Writing "<<str_plot<<" "<<"["<<ymin<<":"<<ymax<<"]"<<std::endl;

    
      stream_plot.open(str_plot);
      //stream_plot<<"set view map\n";
      // stream_plot<<"unset key\n";
      // stream_plot<<"unset tics\n";
      stream_plot<<"unset border\n";
      stream_plot<<"unset colorbox\n";
      stream_plot<<"set yrange ["<<ymin<<":"<<ymax<<"]\n";
      stream_plot<<"set xrange ["<<beginx<<":"<<endx<<"]\n";
      stream_plot<<"set terminal png size 1200, 600\n";
      for (uint i = start; i <= stop; ++i) {
	stream_plot<<"set output \""<<file_name<<i<<".png\"\n";
	stream_plot<<"plot '"<<file_name<<i<<".txt"<<"' with line\n";//, '"<<file_name<<"max_"<<i<<".txt' with line"<<"\n";
	if (acc_max) {
	  stream_plot<<"set output \""<<file_name<<"max_"<<i<<".png\"\n";
	  stream_plot<<"plot '"<<file_name<<"max_"<<i<<".txt"<<"'\n";
	}
      }
      stream_plot<<"set term pop; set out;";
      stream_plot.close();
    }
    if (merge_) {
      std::ofstream stream_plot;
      std::stringstream ss_plot; 
      ss_plot <<file_name<<"_merge"<<".plot";
      std::string str_plot(ss_plot.str());
      //    std::cout<<"Writing "<<str_plot<<std::endl;
    
      stream_plot.open(str_plot);
      //  stream_plot<<"set view map\n";
      //    stream_plot<<"unset key\n";
      //   stream_plot<<"unset tics\n";
      stream_plot<<"unset border\n";
      stream_plot<<"unset colorbox\n";
      stream_plot<<"set yrange ["<<ymin<<":"<<ymax<<"]\n";
      stream_plot<<"set terminal png size 1200, 600\n";
      for (uint i = start; i <= stop; ++i) {
	stream_plot<<"set output \""<<file_name<<"_merge"<<i<<".png\"\n";
	stream_plot<<"plot '"<<file_name<<i<<".txt"<<"' with line, '"<<merge_name<<i<<".txt"<<"' with line\n";
	//      stream_plot<<"replot '"<<merge_name<<i<<".txt"<<"' with line\n";
      }
      stream_plot<<"set term pop; set out;";
      stream_plot.close();
    }
    if (fft_) {
      ip.plotSpectrogram();
      ip.plotSpectrum();
    }
  
  }
  return 0;
}



