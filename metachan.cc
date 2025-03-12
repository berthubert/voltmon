#include <stdio.h>
#include <stdint.h>
#include "nlohmann/json.hpp"
#include <sys/time.h>
#include "misc.hh"

int main()
{
  struct msmt {
    int32_t l;
    int32_t r;
  };
  static_assert(sizeof(msmt)==8);

  nlohmann::json j = nlohmann::json::object();
  j["datarate"]=192000;
  j["channels"]=1;
  j["format"]="S32_LE";
  j["source"]="power";
  j["starttime"]=getSubsecUnixTime();
  printf("%s\n", j.dump().c_str());
  fflush(stdout);
  msmt m;
  while(fread(&m, 1, sizeof(m), stdin) == 8) {
    if(fwrite(&m.r, 1, 4, stdout) != 4)
      break;
  }
    
}
