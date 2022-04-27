#include "definitions.hpp"
#include "FFT.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <vector>

#define WINDOW_FUNC 0

std::string file_name;
int window_size = 1024;
FLOAT dt = 100.0/2048.0;
FLOAT sample_rate = 1/dt;
uint nb_frequencies = 1/(dt*window_size);
FLOAT frequency_step = window_size;

std::list<FLOAT> samples;
  
FLOAT *spectrum_re;
FLOAT *spectrum_im;

std::list<std::vector<FLOAT> > spectrogram;
std::list<std::vector<FLOAT> > spectrogram_im;
std::list<std::vector<FLOAT> > spectrogram_re;

int t = 0;


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

float computeSpectrum() {
  // std::cout<<"compute spectrum "<<spectrogram.size() <<std::endl;
  FLOAT *in_re = new FLOAT[window_size];
  FLOAT *in_im = new FLOAT[window_size];

  std::list<FLOAT>::iterator it;
  uint i = 0;
  for (it = samples.begin(); /*it != samples.end()*/ i < window_size; ++it, ++i) {
    in_re[i] = (*it);
    in_im[i] = 0.0;
  }
  //    WindowFunc(WINDOW_FUNC, window_size, in_re);
  FFT(window_size, false, in_re, in_im, spectrum_re, spectrum_im);
  float max = 0;
  int imax = 0;
  for (int k = 0; k < window_size/2; ++k) {
      //      std::cout<< k*frequency_step<<" "<<k<<std::endl;
     if ((float)k*frequency_step < 0.025) {
       spectrum_re[k] = 0;
       spectrum_im[k] = 0;
       spectrum_re[window_size-1-k] = 0;
       spectrum_im[window_size-1-k] = 0;
     }
    float e = pow(spectrum_re[k], 2) + pow(spectrum_im[k], 2);
    if (e > max) {
      max = e;
      imax = k;
    }

  }
  std::cout<<"max "<<imax*frequency_step<<" "<<imax<<" "<<max<<std::endl;
  // for (i = 0; i < window_size; ++i) {
  //   in_re[i] = 0.0;
  //   in_im[i] = 0.0;
  // }
  //std::cout<<"inverse FFT "<<std::endl;
  FFT(window_size, true, spectrum_re, spectrum_im, in_re, in_im);
  samples.clear();
  //std::cout<<" win size "<<window_size<<std::endl;
  
  for (int k = 0; k < window_size; ++k) {
    samples.push_back(in_re[k]);
    //    samples.pop_back();
  }

  //RealFFT(window_size, in_re, spectrum_re, spectrum_im);
  //FFT(window_size, true, spectrum_re, spectrum_im, in_re, in_im);
  
  std::vector<FLOAT> power_spectrum(nb_frequencies);
  std::vector<FLOAT> spec_re(nb_frequencies);
  std::vector<FLOAT> spec_im(nb_frequencies);
  FLOAT epsilon_db = 1e-18;
  for (uint i = 0; i < nb_frequencies; ++i) {
    FLOAT e = powf(spectrum_re[i]/window_size, 2) + powf(spectrum_im[i]/window_size, 2);
    
    FLOAT e_db = - 10*log10(e+epsilon_db);
    power_spectrum[i] = e;
    spec_re[i] = spectrum_re[i]/window_size;
    spec_im[i] = spectrum_im[i]/window_size;
    // std::cout<<"s "<<i<<" "<<spectrum_re[i]<<" "<<spectrum_im[i]<<std::endl;
    // spectrum_re[i] = in_re[i];
    // spectrum_im[i] = 0;
  }
  // spectrogram.push_back(power_spectrum);
  // spectrogram_re.push_back(spec_re);
  // spectrogram_im.push_back(spec_im);
  // if (spectrogram.size() > 256) {
  //   // plotSpectrum();
  //   // plotSpectrogram();
  //   // plotSamples();
  //   //    exit(0);
  //   spectrogram.pop_front();
  // }
  
  delete[] in_re;
  delete[] in_im;
  return imax*frequency_step;
}


COMPLEX amplitude(FLOAT frequency) {
  uint index = frequency/frequency_step; //TODO interpolation ?
  COMPLEX out(spectrum_re[index], spectrum_im[index]);
  return out;
}

void plotSpectrum(int ind) {
  std::stringstream ss;
  ss <<file_name<<"spectrum_"<<ind<<".txt";
  std::string str(ss.str());
  //    std::cout<<"Plotting "<<str<<std::endl;
  std::ofstream  out_file;
  out_file.open(str.c_str());

  FLOAT epsilon_db = 1e-18; 
  
  for (uint i = 1; i < nb_frequencies/2; ++i) {
    FLOAT e = powf(spectrum_re[i]/window_size, 2) + powf(spectrum_im[i]/window_size, 2);    
    FLOAT e_db = - 10*log10(e  + epsilon_db);
    out_file << i*frequency_step << " " << e <<"\n";
  }
  out_file.close();

  // std::stringstream ss2;
  // ss2 <<file_name<<ind<<"_spectrum_imag.txt";
  // std::string str2(ss2.str());
  // // std::cout<<str<<std::endl;
  // std::ofstream  out_file2;
  // out_file2.open(str2.c_str());
  
  // for (uint i = 1; i < nb_frequencies/2; ++i) {
  //   out_file2 << i*frequency_step << " " << spectrum_re[i]/window_size <<" "<<spectrum_im[i]/window_size <<"\n";
  // }
  // out_file2.close();

    
}

void plotSpectrogram() {
  std::stringstream ss;
  ss <<file_name<<"_spectrogram.txt";
  std::string str(ss.str());
  std::cout<<"Plotting "<<str<<" "<<dt<<std::endl;
  std::ofstream  out_file;
  out_file.open(str.c_str());

  std::list<std::vector<FLOAT> >::const_iterator it;
  std::list<std::vector<FLOAT> >::const_iterator it_re = spectrogram_re.begin();
  std::list<std::vector<FLOAT> >::const_iterator it_im = spectrogram_im.begin();
  uint i = 0;
  for (it = spectrogram.begin(); it != spectrogram.end(); ++it, ++i, ++it_re, ++it_im) {
    for (uint j = 4; j < nb_frequencies/2; ++j) {
      float f_cur = j*frequency_step;
      float k_cur = pow((2*M_PI*f_cur), 2)/9.81;
      float w_cur = 2*M_PI/k_cur;
      out_file << (float)i*dt<<" "<<j*frequency_step<<" "<<w_cur << " " << (*it)[j]/* <<" "<<(*it_re)[j]<<" "<<(*it_im)[j]*/<<"\n";
      // std::cout<<(float)i*dt<<" "<<j*frequency_step<<" "<<w_cur << " " << (*it)[j]/* <<" "<<(*it_re)[j]<<" "<<(*it_im)[j]*/<<"\n";
    }
    out_file <<"\n";
  }
  //    std::cout<<" "<<std::endl;
  out_file.close();
}
void plotSamples(int ind) {
  std::stringstream ss;
  ss <<file_name<<"corr_"<<ind<<".txt";
  std::string str(ss.str());
  std::ofstream  out_file;
  out_file.open(str.c_str());
  
  std::list<FLOAT>::const_iterator it;
  uint i = 0;
  for (it = samples.begin(); it != samples.end(); ++it, ++i) {
    out_file << i*dt<<" "<<(*it) <<"\n";
  }

  out_file.close();
}




int main(int argc, char **argv) {

  int stop = 0, start = 0;
  float endx = 0, beginx = 0;
  float ymin = 0, ymax = 0;
  float ymin_s = 0, ymax_s = 0;
  bool plot_ = false;
  float max_freq = 0;
  int nmax = 0;
  
  //// PARSING  /////
  for (int i = 1;  i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "-file") {
      assert(argc > i+1);
      file_name = argv[i+1];
      std::cout<<"filename: "<<file_name<<std::endl;
      ++i;
    }
        if (s == "-plot") {
      std::cout<<"PLOT"<<std::endl;
      plot_ = true;
    }
    if (s == "-w") {
      assert(argc > i+1);
      window_size = atoi(argv[i+1]);
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
  }
  /////
  
  sample_rate = 1.0/dt;
  frequency_step = 1/(dt*window_size);
  nb_frequencies = window_size;

  // std::cout<<"freq step "<<frequency_step<<std::endl;
  // std::cout<<"nb freq "<<window_size<<std::endl;
  // std::cout<<"freq range "<<window_size*frequency_step<<std::endl;
  // std::cout<<"time range "<<window_size*dt<<std::endl;
  
  for (int ind  = start; ind < stop; ++ind) {
    
    samples.clear();
    for (uint i  = 0; i < window_size; ++i) {
      samples.push_back(0.0);
    }
    spectrum_re = new FLOAT[nb_frequencies];
    spectrum_im = new FLOAT[nb_frequencies];
    for (uint i = 0; i < nb_frequencies; ++i) {
      spectrum_re[i] = 0;
      spectrum_im[i] = 0;
    }

    std::stringstream ss; 
    ss <<file_name<<ind<<"_samples.txt";
    std::string str(ss.str());
    std::ifstream  file;
    file.open(str);
    
    if (!file.good()) {
      std::cerr<<"FILE NOT GOOD: "<<file_name<<std::endl;
      file.close();
      //      return -1;
    } else { 
      std::cout<<"FILE OPEN: "<<file_name<<std::endl;
      int i = 0;
      while (!file.eof()) {
	float x, h;
	file >> x >> h;
	if (endx == 0) {
	  endx = x;
	}
	if (h > ymax) {
	  ymax = h;
	}
	if (h < ymin) {
	  ymin = h;
	}
	//      h = exp(-0.00001*((float)i*dt)*((float)i*dt))*cos(0.6*(float)i*dt*2.0f*(float)M_PI);
	//      std::cout<<h<<" ";
	samples.push_front(h);
	samples.pop_back();
	if (t == 2000) {
	  //	  std::cout<<"samples size "<<samples.size()<<std::endl;
	  float m = computeSpectrum();
	  max_freq += m;
	  ++nmax;
	  for (uint i = 1; i < nb_frequencies/2; ++i) {
	     FLOAT e = powf(spectrum_re[i]/window_size, 2) + powf(spectrum_im[i]/window_size, 2);    
	     //  FLOAT e_db = - 10*log10(e  + epsilon_db);
	     if (e > ymax_s) {
	       ymax_s = e;
	     }
	     if (e < ymin_s) {
	       ymin_s = e;
	     }
	  }
	  //  std::cout<<"samples size --- "<<samples.size()<<" "<<ind<<std::endl;
	  t = 0;
	  break;
	}
	++i;
	++t;
      }
    }
    plotSamples(ind);
    plotSpectrum(ind);
    //  plotSpectrogram();
    clear();
  }
  if (plot_) {
    std::ofstream stream_plot;
    std::stringstream ss_plot; 
    ss_plot <<file_name<<"corr.plot";
    std::string str_plot(ss_plot.str());
    std::cout<<"Writing "<<str_plot<<" "<<"["<<ymin<<":"<<ymax<<"]"<<std::endl;
    
    stream_plot.open(str_plot);
    if (!stream_plot.good()) {
      std::cerr<<str_plot<<" cannot be opened."<<std::endl;
      return -1;
    }
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
      stream_plot<<"plot '"<<file_name<<"corr_"<<i<<".txt"<<"' with line\n";//, '"<<file_name<<"max_"<<i<<".txt' with line"<<"\n";
      // if (acc_max) {
      // 	stream_plot<<"set output \""<<file_name<<"max_"<<i<<".png\"\n";
      // 	stream_plot<<"plot '"<<file_name<<"max_"<<i<<".txt"<<"'\n";
      // }
    }
    stream_plot<<"set term pop; set out;";
    stream_plot.close();

    float beginx_s = 0, endx_s = 1.5;//window_size/2*frequency_step;
        
    std::ofstream stream_plot_s;
    std::stringstream ss_plot_s;
    ss_plot_s <<file_name<<"spec.plot";
    std::string str_plot_s(ss_plot_s.str());
    std::cout<<"Writing "<<str_plot_s<<" "<<"["<<ymin<<":"<<ymax<<"]"<<std::endl;

    
    stream_plot_s.open(str_plot_s);
    //stream_plot<<"set view map\n";
    // stream_plot<<"unset key\n";
    // stream_plot<<"unset tics\n";
    stream_plot_s<<"unset border\n";
    stream_plot_s<<"unset colorbox\n";
    stream_plot_s<<"set yrange ["<<ymin_s<<":"<<ymax_s<<"]\n";
    stream_plot_s<<"set xrange ["<<beginx_s<<":"<<endx_s<<"]\n";
    stream_plot_s<<"set terminal png size 1200, 600\n";
    for (uint i = start; i <= stop; ++i) {
      stream_plot_s<<"set output \""<<file_name<<"spec_"<<i<<".png\"\n";
      stream_plot_s<<"plot '"<<file_name<<"spectrum_"<<i<<".txt"<<"' with line\n";//, '"<<file_name<<"max_"<<i<<".txt' with line"<<"\n";
      // if (acc_max) {
      // 	stream_plot<<"set output \""<<file_name<<"max_"<<i<<".png\"\n";
      // 	stream_plot<<"plot '"<<file_name<<"max_"<<i<<".txt"<<"'\n";
      // }
    }
    stream_plot_s<<"set term pop; set out;";
    stream_plot_s.close();
  }
  std::cout<<"Mean max "<<max_freq/(float)nmax<<std::endl;
  return 0;
}
