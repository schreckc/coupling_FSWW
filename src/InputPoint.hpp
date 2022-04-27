#ifndef INPUTPOINT_HPP
#define INPUTPOINT_HPP

#include <list>
#include <vector>
#include "definitions.hpp"
#include "FFT.hpp"
#include <fstream>
#include <iostream>

/*
 * 0: Rectangular (no window)
 * 1: Bartlett    (triangular)
 * 2: Hamming
 * 3: Hanning
 * 4: Blackman
 * 5: Blackman-Harris
 * 6: Welch
 * 7: Gaussian(a=2.5)
 * 8: Gaussian(a=3.5)
 * 9: Gaussian(a=4.5)
 */
#define WINDOW_FUNC 0

class InputPoint {
private:
  VEC2 pos;
  
  uint window_size;
  FLOAT dt;
  FLOAT sample_rate;
  uint nb_frequencies;
  FLOAT frequency_step;
  
  std::list<FLOAT> samples;
  std::string name;

  float ymax;

  
public:
  std::vector<float> maxf, maxf2, meanf, middlef, ampli_maxf, ampli_maxf2, total_ampli;

  
  FLOAT *spectrum_re;
  FLOAT *spectrum_im;

  std::list<std::vector<FLOAT> > spectrogram;
  std::list<std::vector<FLOAT> > spectrogram_im;
  std::list<std::vector<FLOAT> > spectrogram_re;

  uint t;

  bool rec;
    
public:
  InputPoint() {
    pos = VEC2(0, 0);
    window_size = 0;
    dt = 0;
    nb_frequencies = 0;
    frequency_step = 0;
    
    spectrum_re = NULL;
    spectrum_im = NULL;
    
    t = 0;
    name = "test";
    rec = false;
    ymax = 0;
    //std::cout<<"Win size 000"<<window_size<<" "<<samples.size()<<std::endl;
  }
  InputPoint(uint nb_samples, FLOAT time_step, bool r = false) {
    pos = VEC2(0, 0);
    window_size = nb_samples;
    dt = time_step;
    sample_rate = 1.0/dt;
    frequency_step = 1/(dt*window_size);
    nb_frequencies = window_size/2;
    samples.clear();
    
    std::cout<<"freq step "<<frequency_step<<std::endl;
    std::cout<<"nb freq "<<window_size<<std::endl;
    std::cout<<"freq range "<<window_size*frequency_step<<std::endl;
    std::cout<<"time range "<<window_size*dt<<std::endl;

    rec = r;
    
    if (!rec) {
      std::cout<<"REC NOT"<<std::endl;
     for (uint i  = 0; i < window_size; ++i) {
       samples.push_back(0.0);
     }
    }
     if (rec) {
       std::cout<<"REC "<<samples.size()<<std::endl;
     }

    
     std::cout<<"Win size "<<window_size<<" "<<samples.size()<<std::endl;
     spectrum_re = new FLOAT[/*nb_frequencies*/window_size];
    spectrum_im = new FLOAT[/*nb_frequencies*/window_size];
    for (uint i = 0; i < /*nb_frequencies*/window_size; ++i) {
      spectrum_re[i] = 0;
      spectrum_im[i] = 0;
    }
    t = 0;
    name = "aaarg";
    ymax = 0;
  }
    
  ~InputPoint() {
     // if (spectrum_re != NULL) {
     //   delete[]spectrum_re;
     //   spectrum_re = NULL;
     // }
     // if (spectrum_im != NULL) {
     //   delete[]spectrum_im;
     //   spectrum_im = NULL;
     // }
    // plotSpectrum();
    // plotSpectrogram();
    // plotSamples();
  }

  void clear() {
   if (spectrum_re != NULL) {
       delete[]spectrum_re;
       spectrum_re = NULL;
     }
     if (spectrum_im != NULL) {
       delete[]spectrum_im;
       spectrum_im = NULL;
     }
  }
  
  void setPos(VEC2 p) {
    pos = p;
  }
  void setPos(FLOAT x, FLOAT y) {
    pos = VEC2(x, y);
  }

  VEC2 getPos() const {
    return pos;
  }

  void setName(std::string n) {
    name = n;
  }
  
  void update(FLOAT next_sample) {
        
     samples.push_front(next_sample);
     if (!rec && samples.size() > window_size) {
       //  std::cout<<"REC NOT"<<std::endl;
       samples.pop_back();
     }
     //      if ((t + 1) % (window_size) == 0 /*|| t == 1999*/) {
      //     if (t%32 == 31) {
      //     if (t > 1) {
        computeSpectrum();
        //   t = 0;
	//  }
      ++t;
      //std::cout<<"NEXT sample "<<next_sample<<"  "<<t<<"   "<<samples.size()<<" "<<nb_frequencies<<std::endl;   
    
  }
  
  void computeSpectrum()  {
    //std::cerr<<"compute spectrum "<<spectrogram.size()<<" "<<t <<" "<<nb_frequencies<<std::endl;
    if (nb_frequencies > 0) {
    FLOAT *in_re = new FLOAT[window_size];
    FLOAT *in_im = new FLOAT[window_size];

    std::list<FLOAT>::iterator it;
    uint i = 0;
    for (it = samples.begin(); /*it != samples.end()*/ i < window_size;  ++i) {
      if (it != samples.end()) {
	in_re[i] = (*it);
	++it;
      } else {
	in_re[i] = 0.0;
      }
      in_im[i] = 0.0;
    }
    //   std::cout<<"WIN size IP  -- "<<window_size<<std::endl;
    //    WindowFunc(WINDOW_FUNC, window_size, in_re);
    for (uint i = 0; i < /*nb_frequencies*/window_size; ++i) {
      spectrum_re[i] = 0;
      spectrum_im[i] = 0;
    }

    
    FFT(window_size, false, in_re, in_im, spectrum_re, spectrum_im);
    //RealFFT(window_size, in_re, spectrum_re, spectrum_im);
    //FFT(window_size, true, spectrum_re, spectrum_im, in_re, in_im);
  
    std::vector<FLOAT> power_spectrum(nb_frequencies);
    std::vector<FLOAT> spec_re(nb_frequencies);
    std::vector<FLOAT> spec_im(nb_frequencies);
    FLOAT epsilon_db = 1e-18;
    float max = 0;
    int imax = 0;
    float max2 = 0;
    int imax2 = 0;
    float total_a = 0;
    float mean = 0;
    for (uint i = 3; i < nb_frequencies; ++i) {
      // WARNING added used win_size/2 instead of win_size
      FLOAT e = powf(spectrum_re[i]/(float)window_size*2.0, 2) +
	powf(spectrum_im[i]/(float)window_size*2.0, 2);
    
      FLOAT e_db = - 10*log10(e+epsilon_db);
      power_spectrum[i] = e;
      if (e > ymax) {
	ymax = e;
      }
      spec_re[i] = spectrum_re[i]/(window_size/2.0);
      spec_im[i] = spectrum_im[i]/(window_size/2.0);
      //      std::cout<<"s "<<i<<" "<<spectrum_re[i]<<" "<<spectrum_im[i]<<" "<<power_spectrum[i]<<std::endl;
      // spectrum_re[i] = in_re[i];
      // spectrum_im[i] = 0;
      total_a += e;
      //      std::cout<<"ampli "<<total_a<<" "<<e<<std::endl;
      mean += e*i;
      if (e > max && i > 0) {
      	max2 = max;
      	imax2 = imax;
      	max = e;
      	imax = i;
      } else if (e > max2 && i > 0) {
      	max2 = e;
      	imax2 = i;
      }
    }
    //    std::cout<<"YMAX "<<ymax<<" "<<t<<std::endl;
     // std::cout<<"\nWinsize: "<< window_size<<"  frame: "<<t<<std::endl;
     // std::cout<<"AMPLI: "<< total_a<<std::endl;
     // std::cout<<"MAX freq "<<imax*frequency_step<<" "<<imax<<" "<<max<<" "<<max/total_a<<" "<<total_a<<std::endl;
     //  std::cout<<"Mean freq "<< mean/total_a*frequency_step<<std::endl;
     //  std::cout<<"Middle freq "<< (max*imax+max2*imax2)/(max+max2)*frequency_step<<std::endl;
     //  std::cout<<"MAX freq 2 "<<imax2*frequency_step<<" "<<imax2<<" "<<max2<<" "<<max2/total_a<<" "<<total_a<<std::endl;
     //  std::cout<<"\n"<<std::endl;

     maxf.push_back(imax*frequency_step);
     maxf2.push_back(imax2*frequency_step);
     meanf.push_back(mean/total_a*frequency_step);
     middlef.push_back((max*imax+max2*imax2)/(max+max2)*frequency_step);
     ampli_maxf.push_back(max);
     ampli_maxf2.push_back(max2);
     total_ampli.push_back(total_a);
     
    if (t%32 == 31) {
      spectrogram.push_back(power_spectrum);
      spectrogram_re.push_back(spec_re);
      spectrogram_im.push_back(spec_im);
      if (!rec && spectrogram.size() > 256) {
    	//
    	// plotSpectrogram();
    	// plotSamples();
    	//    exit(0);
    	spectrogram.pop_front();
      }
    }
    //     plotSamples();
    // if (t%256 == 255) {
    //   plotSpectrum();
    // }
    delete[] in_re;
    delete[] in_im;
    }
  }


  COMPLEX amplitude(FLOAT frequency) const  {
    uint index = frequency/frequency_step; //TODO interpolation ?
    COMPLEX out(spectrum_re[index], spectrum_im[index]);
    return out;
  }

  void plotSpectrum() const {
    std::stringstream ss;
    ss <<name<<"_spectrum/sp"<<t<<".txt";
    std::string str(ss.str());
    std::cout<<"Plotting "<<str<<std::endl;
    std::ofstream  out_file;
    out_file.open(str.c_str());
    float amp = 0;

     std::cout<<"Plotting "<<str<<std::endl;
    
    FLOAT epsilon_db = 1e-18; 
  
    for (uint i = 1; i < nb_frequencies; ++i) {
      FLOAT e = powf(spectrum_re[i]/((float)window_size/2), 2) + powf(spectrum_im[i]/((float)window_size), 2);
      e = sqrt(e);
      amp += e;
      FLOAT e_db = - 10*log10(e  + epsilon_db);
      out_file << i*frequency_step << " " << e <<"\n";
    }
    out_file.close();

    std::cout<< "AMPLI plot "<<amp<<std::endl;
    
    // std::stringstream ss2;
    // ss2 <<name<<"_spectrum_imag.txt";
    // std::string str2(ss2.str());
    // // std::cout<<str<<std::endl;
    // std::ofstream  out_file2;
    // out_file2.open(str2.c_str());
  
    // for (uint i = 1; i < nb_frequencies/2; ++i) {
    //   out_file2 << i*frequency_step << " " << spectrum_re[i]/window_size <<" "<<spectrum_im[i]/window_size <<"\n";
    // }
    // out_file2.close();

    
  }

  void plotSpectrogram() const {
    std::stringstream ss;
    ss <<name<<"_spectrogram.txt";
    std::string str(ss.str());
    std::cout<<"Plotting "<<str<<" "<<dt<<" "<<spectrogram.size()<<std::endl;
    std::ofstream  out_file;
    out_file.open(str.c_str());

    std::list<std::vector<FLOAT> >::const_iterator it;
    std::list<std::vector<FLOAT> >::const_iterator it_re = spectrogram_re.begin();
    std::list<std::vector<FLOAT> >::const_iterator it_im = spectrogram_im.begin();
    uint i = 0;
    uint n = spectrogram.size();
    for (it = spectrogram.begin(); it != spectrogram.end(); ++it, ++i, ++it_re, ++it_im) {
      for (uint j = 1; j < nb_frequencies; ++j) {
	float f_cur = j*frequency_step;
	float k_cur = pow((2*M_PI*f_cur), 2)/9.81;
	float w_cur = 2*M_PI/k_cur;
	out_file << (float)i*dt*32<<" "<<j*frequency_step<<" "<<w_cur << " " << (*it)[j]/* <<" "<<(*it_re)[j]<<" "<<(*it_im)[j]*/<<"\n";
	// std::cout<<(float)i*dt<<" "<<j*frequency_step<<" "<<w_cur << " " << (*it)[j]/* <<" "<<(*it_re)[j]<<" "<<(*it_im)[j]*/<<"\n";
      }
      out_file <<"\n";
    }
    //    std::cout<<" "<<std::endl;
    out_file.close();
  }
  void plotSamples() const {
    std::stringstream ss;
    ss <<name<<"_samples.txt";
    std::string str(ss.str());
    std::ofstream  out_file;
    out_file.open(str.c_str());
  
    std::list<FLOAT>::const_iterator it;
    uint i = 0;
    uint n = samples.size();
    //     std::cout<<"sample size "<<samples.size()<<std::endl;
    for (it = samples.begin(); it != samples.end(); ++it, ++i) {
      out_file << (n-i-1)*dt<<" "<<(*it) <<"\n";
    }

    out_file.close();

    std::stringstream ss_p;
    ss_p <<name<<"_spectrum"<<".plot";
    std::string str_p(ss_p.str());
    std::ofstream  plot_file;
    plot_file.open(str_p.c_str());
    plot_file<<"unset border\n";
    plot_file<<"unset colorbox\n";
    plot_file<<"set yrange ["<<0<<":"<<ymax<<"]\n";
    plot_file<<"set terminal png size 1200, 600\n";
    for (uint i = 1; i < t; ++i) {
      plot_file << "set output \""<<name<<"_spectrum/sp"<<i<<".png\"\n";
      plot_file<<"plot '"<<name<<"_spectrum/sp"<<i<<".txt' with line\n";
    }
    plot_file<<"set term pop; set out;";
    plot_file.close();
  }

  void move(VEC2 trans) {
    pos += trans;
  }

  std::ofstream & record_samples(std::ofstream & file) {
    file << samples.back() <<" ";
    return file;
  }

};



#endif




  // //FFT
  // #include "FFT.h"
  //   float *out = new float[size_slide];
  //   float *out2 = new float[size_slide];
  //   float *in = new float[size_slide];
  //   float *in2 = new float[size_slide];

  //    for (int i = 0; i < size_slide; ++i) {
  //     float f = i * freq_step;
  //     double s = /*(1 + 1*sp.speed * exp(0.001*i)) */ sqrt(pow(r_spectrum[i], 2) + pow(i_spectrum[i], 2));
  //     double phase = 2*M_PI*(rand() / (RAND_MAX + 1.0) - 0.5);
  //     in[i] = cos(phase)*s;
  //     in2[i] = sin(phase)*s;
  //     // if (in[i] != 0 || in2[i] != 0) {
  //     //   INFO("freq reel im", f<<" "<<in[i]<<" "<<in2[i]);
  //     // }
  //   }
  //   FFT(size_slide, true, in, in2, out, out2);

