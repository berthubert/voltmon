#include <stdio.h>
#include <stdint.h>

int main()
{
  struct msmt {
    int32_t l;
    int32_t r;
  };
  static_assert(sizeof(msmt)==8);
  msmt m;
  while(fread(&m, 1, sizeof(m), stdin) == 8) {
    if(fwrite(&m.r, 1, 4, stdout) != 4)
      break;
  }
    
}
