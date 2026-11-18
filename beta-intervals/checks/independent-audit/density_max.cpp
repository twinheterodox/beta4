#include <algorithm>
#include <array>
#include <cstdint>

static int64_t down(int64_t n, int64_t d) {return n/d-(n%d<0);}
static int64_t residue(int64_t n, int64_t d) {return (n%d+d)%d;}
extern "C" int envelope(int n, int b) {
    int64_t d=96LL*b, maxy= n<2*b ? 2*96LL*(n-b) : 2*d;
    std::array<int64_t,66> points{};
    int h=0;
    points[h++]=0;
    points[h++]=d;
    for (int j=0;j<=30;++j){
        points[h++]=residue((int64_t)j*n,d);
        points[h++]=residue((96LL-j)*n,d);
    }
    points[h++]=residue(-72LL*n-d/2,d);
    points[h++]=residue(168LL*n-d/2,d);
    // Both the value at each exact breakpoint and immediately to its right.
    int best=0;
    for (int64_t corner : points) for (int side=0;side<=1;++side){
        int64_t y=2*corner+side;
        if (y>maxy) continue;
        int f=2-3*(down(y+144LL*n+d,2*d)-down(y-336LL*n+d,2*d));
        for (int j=0;j<=30;++j)
            f+=down((192LL-2*j)*n-y,2*d)+down(y-2LL*j*n,2*d)+1;
        best=std::max(best,f);
    }
    return 16*best-221*(n/b);
}
