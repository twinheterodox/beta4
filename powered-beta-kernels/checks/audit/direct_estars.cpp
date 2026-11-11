#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include <cstdint>
using I = int64_t;
I floor_div(I a,I b){I q=a/b,r=a%b;return q-(r<0);}
I mod(I a,I b){I r=a%b;return r<0?r+b:r;}
int main(){
 std::vector<std::pair<I,I>> blocks{{0,8}};
 for(I k=1;k<=8;k++)blocks.push_back({3*k,1});
 for(I k=10;k<=46;k+=2)blocks.push_back({3*k,1});
 I a,b;
 while(std::cin>>a>>b){
  I den=480*b, common=2*den;
  std::vector<I> cuts;
  for(auto [s,m]:blocks){cuts.push_back(mod(s*a,den));cuts.push_back(mod((480-s)*a,den));}
  cuts.push_back(mod(-560*a-240*b,den));
  cuts.push_back(mod(1040*a-240*b,den));
  std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
  I best=INT64_MIN;
  auto evaluate=[&](I y){
   I value=2;
   for(auto [s,m]:blocks){
    // Closed interval count: floor(high-y)-ceil(low-y)+1.
    value+=m*(floor_div(2*(480-s)*a-y,common)+floor_div(y-2*s*a,common)+1);
   }
   // Odd integers in (2y-2(1+r)x,2y+2rx].
   value-=3*(floor_div(y+1120*a+480*b,common)-floor_div(y-2080*a+480*b,common));
   best=std::max(best,value);
  };
  for(size_t j=0;j<cuts.size();j++){
   evaluate(2*cuts[j]);
   I next=j+1<cuts.size()?cuts[j+1]:cuts[0]+den;
   evaluate(cuts[j]+next);
  }
  std::cout<<10*best-179*(a/b)<<'\n';
 }
}
