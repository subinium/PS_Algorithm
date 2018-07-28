typedef unsigned long long T;

T modbase;

T mul(T a, T b){
    if(b==0) return 0;
    T ret=mul(a,b>>1);
    ret=(ret+ret)%modbase;
    if(b&1) ret=(ret+a)%modbase;
    return ret;
}

T pow(T a, T b){
    if(b==0) return 1;
    T ret=pow(a,b>>1);
    ret=mul(ret,ret);
    if(b&1) ret=mul(ret,a);
    return ret;
}

bool miller_counter(T x,T base)
{
    if(base%x==0) return 0;
    int i;
    T d=x-1;
    for(;;){
        if(d&1){
            T tmp=pow(base,d);
            if(tmp!=1 && tmp!=x-1) return 1;
            else return 0;
        } else {
            if(pow(base,d)==x-1) return 0;
            d>>=1;
        }
    }
}

int miller_base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,0};

bool isprime(T x) {
    modbase=x;
    int i;
    for(i=0; i<12; ++i) if(miller_counter(x,miller_base[i])) return 0;
    return 1;
}
