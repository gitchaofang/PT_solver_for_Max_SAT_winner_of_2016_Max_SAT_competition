#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
long seedgen(void){
  long s, seed, pid;
  time_t seconds;
  pid = getpid();
  s = time ( &seconds ); /* get CPU seconds since 01/01/1970 */
  seed = abs(((s*181)*((pid-83)*359))%104729);
  return seed;
}
