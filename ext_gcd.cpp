using namespace std;


// ax + by = gcd(a,b) => x ? y ?
pair<int,int> ext_gcd(int a,int b){
    if(b){
        auto tmp = ext_gcd(b, a%b);
        return {tmp.second, tmp.first - (a/b) * tmp.second};
    }
    else return {1, 0};
}


// ax = 1 mod M , x?
int mod_inv(int a, int M){
    return (ext_gcd(a, M).first + M) % M;
}
