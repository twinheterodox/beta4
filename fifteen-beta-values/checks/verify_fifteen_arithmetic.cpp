// A second implementation: evaluate every floor count directly at the
// cells in y, rather than using the first checker's cumulative jump sweep.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using I=std::int64_t;
I floordiv(I a,I b){I q=a/b,r=a%b;return q-(r<0);}
I mod(I a,I b){I r=a%b;return r<0?r+b:r;}
int main(int argc,char**argv){
 if(argc!=2)return 2;
 std::ifstream f(argv[1]);if(!f)throw std::runtime_error("missing trace");
 I a,b,expected_num,expected_den;std::size_t checked=0;
 while(f>>a>>b>>expected_num>>expected_den){
  I den=96*b,cap=a<2*b?96*(a-b):den;
  std::vector<I> pts{0,cap};
  for(int j=0;j<=30;j++){
   for(I z:{j*a,(96-j)*a}){I y=mod(z,den);if(y>0&&y<cap)pts.push_back(y);}
  }
  for(I z:{-72*a-48*b,168*a-48*b}){I y=mod(z,den);if(y>0&&y<cap)pts.push_back(y);}
  std::sort(pts.begin(),pts.end());pts.erase(std::unique(pts.begin(),pts.end()),pts.end());
  I best=INT64_MIN;
  for(std::size_t t=1;t<pts.size();t++){
   I y2=pts[t-1]+pts[t],D=2*den,v=0;
   for(int j=0;j<=30;j++)v+=floordiv(y2-2*j*a,D)-floordiv(y2-2*(96-j)*a,D);
   v-=3*(floordiv(y2+144*a+96*b,D)-floordiv(y2-336*a+96*b,D));
   best=std::max(best,v);
  }
  I sixteen_e=16*std::max<I>(0,best+2)-221*(a/b);
  if(sixteen_e*expected_den!=16*expected_num){
   std::cerr<<"Mismatch "<<a<<"/"<<b<<" got "<<sixteen_e<<"/16 expected "<<expected_num<<"/"<<expected_den<<"\n";
   return 1;
  }
  checked++;
 }
 std::cout<<"{\"status\":\"passed\",\"direct_floor_intervals\":"<<checked<<"}\n";
}
