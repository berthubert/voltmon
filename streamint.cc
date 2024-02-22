#include "fmt/core.h"
#include "fmt/os.h"
#include <stdexcept>
#include <math.h>
#include <iostream>
#include <vector>
#include <time.h>
#include "sqlwriter.hh"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>

using namespace boost::accumulators;
using namespace std;
#include "nr3.h"
#include "amoeba.h"
#include "brent.hh"


/*
Receiver of a stream of best data we can get (32 bit, 192kHz for example).
We do anomaly detection at this prime data rate, since it allows us to see the most.
We have a ringbuffer of the full data

We study every second independently
If there is an anomaly, emit past few seconds at high resolution, and continue emitting for a few seconds more. If there are subsequent anomalies, extend the period.

*/

int32_t getVal(const vector<int32_t>& vals, int hz, double t, double w)
{
  int64_t begin = (t - w/2) * hz;
  int64_t end = (t + w/2) * hz;
  if(begin < 0)
    return *vals.begin();
  
  if(end + 1 >= (int64_t)vals.size())
    return *vals.rbegin();

  if(begin == end)
    return vals.at(begin);
  
  int64_t s = 0;
  for(size_t pos = begin; pos < (size_t)end; ++pos) // this is the dumb average
    s += vals.at(pos);
  return s / (end - begin);
}

vector<int32_t> getVals(const vector<int32_t>& vals, int hz, double begin, double end, double step)
{
  vector<int32_t> ret;
  ret.reserve((end-begin)/step);
  for(double t = begin ; t < end; t+= step)
    ret.push_back(getVal(vals,hz,t, step));
  return ret;
}

double vdiff(const vector<int32_t>& a, const vector<int32_t>& b)
{
  if(a.size() != b.size())
    throw std::runtime_error("Unequal size");
  double diff=0;
  for(size_t pos = 0; pos < a.size(); ++pos)
    diff += abs(b[pos] - a[pos]);
  return diff;
}
double sdiff(const vector<int32_t>& a, const vector<int32_t>& b)
{
  if(a.size() != b.size())
    throw std::runtime_error("Unequal size");
  double diff=0;
  for(size_t pos = 0; pos < a.size(); ++pos)
    diff += pow(b[pos] - a[pos], 2);
  return diff;
}

double sexcess(const vector<int32_t>& a, const vector<int32_t>& b, double lim)
{
  if(a.size() != b.size())
    throw std::runtime_error("Unequal size");

  uint32_t excess = 0;
  for(size_t pos = 0; pos < a.size(); ++pos)
    if(abs(b[pos] - a[pos]) > fabs(lim))
      excess++;
  return 1.0*excess/a.size();
}

double maxdiff(const vector<int32_t>& a, const vector<int32_t>& b)
{
  if(a.size() != b.size())
    throw std::runtime_error("Unequal size");

  uint32_t mdiff=0;
  for(size_t pos = 0; pos < a.size(); ++pos) {
    uint32_t diff = abs(b[pos] - a[pos]);
    if(diff > mdiff)
      mdiff = diff;
  }
  return mdiff;
}


double asum(const vector<int32_t>& a)
{
  double ret = 0;
  for(const auto& v : a)
    ret += abs(v);
  return ret;
}


void getData(FILE* fp, vector<int32_t>& data, uint32_t samples)
{
  data.resize(samples);
  int rc = fread(&data[0], 4, samples, fp);
  if(rc >= 0)
    data.resize(rc);
  else
    throw std::runtime_error("fread: "+string(strerror(errno)));
}

int main(int argc, char** argv)
{
  SQLiteWriter sqlw("voltmon.sqlite");
  fmt::print("t,tstamp,absv,minv,maxv,q01,q99,phz,phz2,relabsdev,dc,amp,offset,percexcess5,percexcess10\n");
  unsigned int hz = 192000;
  vector<int32_t> data;
  typedef accumulator_set<int32_t, stats<tag::p_square_quantile> > bacc_t;
  for(int t=0;;++t) {

    getData(stdin, data, hz);
    if(data.size() != hz) // whole seconds only, discard last bit
      break;

    if(!t) // first second often has bad data, let's skip it
      continue;

    double step=0.0001; // 0.1ms, 10kHz, downsample
    auto interp = getVals(data, hz, 0.1, 0.2, step); // 0.1 seconds for Hz determination

    uint64_t abssum=0;
    int32_t minv= std::numeric_limits<int32_t>::max();
    int32_t maxv= std::numeric_limits<int32_t>::min();
    bacc_t acc0(quantile_probability = 0.01);
    bacc_t acc1(quantile_probability = 0.99);

    for(const auto& v : interp) {
      abssum += abs(v);
      minv = min(v, minv);
      maxv = max(v, maxv);
      acc0(v);
      acc1(v);
    }

    double q01= extract_result< tag::p_square_quantile >( acc0);
    double q99= extract_result< tag::p_square_quantile >( acc1) ;

    
    auto func = [&](double phz) {
      double offset = 1.0/phz;
      auto vals2 = getVals(data, hz, 0.1 + offset, 0.2 + offset, step);
      return fabs(vdiff(interp, vals2));
    };
    Brent b;
    b.ax=35;
    b.bx=50;
    b.cx=65;
    double phz = b.minimize(func);

    interp = getVals(data, hz, 0.0, 1.0, step); // whole chunk
    auto genWave = [&](const vector<double> & params)
    {
      // dc + amp * sin(2*Pi*freq*(t-offset))
      const double& dc = params[0];
      const double& amp = params[1];
      const double& freq = params[2];
      const double& offset = params[3];
      vector<int32_t> gen;
      gen.resize(interp.size());

      //      fmt::print("a {}, b {}\n", a, b);
      //auto out = fmt::output_file("fit-"+to_string(ctr++)+".csv");      
      //out.print("t,gen,interp\n");
      for(size_t pos = 0 ; pos < gen.size(); ++pos) {
        const double t = 1.0*pos*step;
        gen[pos] = dc + amp * sin(2*M_PI*freq*(t-offset));
        //        out.print("{},{},{}\n", t, gen[pos], interp[pos]);
      }
      return gen;
    };
    auto func2= [&](const vector<double>& params)
    {
      auto gen = genWave(params);
      auto ret = vdiff(interp, gen);
      return ret;
    };

    Amoeba amo(0.001);
    double offsethint=0;
    double dc=(maxv+minv)/2.0; // or use q99?
    for(size_t pos = 1 ; pos < 0.1*interp.size(); ++pos) {
      if(interp[pos-1] < dc && interp[pos] > dc) {
        offsethint = pos * step;
        break;
      }
    }
    //                    dc                        amp
    vector<double> params{dc,                   (maxv-minv)/2.0,      phz, offsethint};
    vector<double> deltas{0.1*(maxv-minv)/2.0, 0.1*(maxv-minv)/2.0, 0.01, 0.0001};

    // do the simplex optimization
    vector<double> res = amo.minimize(params, deltas, func2);

    // plot the fit, and see how good it really is
    auto gen = genWave(res);
    double absreldev = sdiff(gen,interp)/asum(gen);
    double percexcess5 = 100.0 * sexcess(gen, interp, 0.05*res[1]); // amplitude, 5%
    double percexcess10 = 100.0 * sexcess(gen, interp, 0.10*res[1]); // amplitude, 10%
    double maxpercdev = 100.0*maxdiff(gen, interp)/res[1];
    time_t now = time(nullptr);
    fmt::print("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n", t, now, 1.0*abssum/data.size(), minv, maxv, q01, q99, phz, res[2],
               absreldev, res[0], res[1], res[3], percexcess5, percexcess10,maxpercdev);

    // fmt::print("t,tstamp,absv,minv,maxv,q01,q99,phz,phz2,relabsdev,dc,amp,offset,percexcess5,percexcess10\n");

    sqlw.addValue({{"t", t}, {"tstamp", now}, {"absv", 1.0*abssum/data.size()}, {"minv", minv}, {"maxv", maxv}, {"q01", q01}, {"q99", q99}, {"phz", phz}, {"phz2", res[2]}, {"absreldev", absreldev}, {"dc", res[0]}, {"amp", res[1]}, {"offset", res[3]}, {"percexcess5", percexcess5}, {"percexcess10", percexcess10}, {"maxpercdev", maxpercdev}});
                  

    fflush(stdout);

    if(absreldev > 370000 || maxv > 1.05 * q99 || minv < 1.05 * q01 || percexcess10 > 1 || percexcess5 > 6 || maxpercdev > 11) {
      fmt::print(stderr, "Anomaly at {}, inmag = {}, dc = {}, amp = {}, phz2 = {}\n", t, params[1], res[0], res[1], res[2]);
      auto out = fmt::output_file("anomaly-"+to_string(t)+".csv");
      out.print("t,gen,orig\n");
      for(size_t pos = 0 ; pos < gen.size(); ++pos) {
        const double t2 = 1.0*pos*step;
        out.print("{},{},{}\n", t2, gen[pos], interp[pos]);
      }
      vector<uint8_t> raw;
      raw.resize(4* data.size());
      raw.assign((uint8_t*)&data[0], ((uint8_t*)&data[0]) + raw.size());
      sqlw.addValue({{"t", t}, {"tstamp", now}, {"hz", hz}, {"format", "S32_LE"}, {"raw", raw}}, "dumps");
    }
  }
}
