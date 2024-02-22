#include "wavint.hh"
#include "fmt/core.h"
#include "fmt/os.h"
#include <stdexcept>
#include <math.h>
#include <iostream>
#include "brent.hh"

using namespace std;

double vdiff(const vector<float>& a, const vector<float>& b)
{
  if(a.size() != b.size())
    throw std::runtime_error("Unequal size");
  double diff=0;
  for(size_t pos = 0; pos < a.size(); ++pos)
    diff += fabs(b[pos] - a[pos]);
  return diff;
}

int main(int argc, char** argv)
{
  WAVInt wi(argv[1]);
  double step = 0.0001; // 10kHz data
  double width = 0.2; // in seconds, so 10 waves
  
  fmt::print("t,f\n");
  cerr<<"Have "<< wi.lengthS() <<" seconds of data"<<endl;

  for(double t=0; t < wi.lengthS(); t += 1.0) {
    auto vals = wi.getVals(t, t + width, step);

    auto func = [&](double hz) {
      double offset = 1.0/hz;
      auto vals2 = wi.getVals(t + offset, t + width + offset, step);
      return fabs(vdiff(vals,vals2));
    };
    Brent b;
    b.ax=35;
    b.bx=50;
    b.cx=65;
    fmt::print("{},{}\n", t, b.minimize(func));
    fflush(stdout);

    /*
    // b.xmin is now set
    double offset = 1/b.xmin; // seconds

    vector<float> mold;
    mold.resize(offset/step);
    int numwaves=10;
    for(unsigned int n = 0 ; n < numwaves; ++n) {
      auto w = wi.getVals(t + n * offset, t + (n+1)*offset, step); 
      for(unsigned int s = 0 ; s < w.size(); ++s)
        mold[s] += vals.at(s);
    }
    for(auto& m : mold)
      m/=numwaves;
    auto out = fmt::output_file(to_string((int)t)+".csv");
    out.print("t,v\n");
    for(unsigned int pos = 0; pos < mold.size(); ++pos)
      out.print("{},{}\n", pos*step, mold.at(pos));
    */
  }

  
}
