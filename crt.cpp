#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

vector<ll> rem;
vector<ll> mod;

// ax + by = gcd(a,b) => x ? y ?
pair<ll,ll> ext_gcd(ll a,ll b){
    if(b){
        auto tmp = ext_gcd(b, a%b);
        return {tmp.second, tmp.first - (a/b) * tmp.second};
    }
    else return {1, 0};
}

// ax = 1 mod M , x?
ll mod_inv(ll a, ll M){
    return (ext_gcd(a, M).first + M) % M;
}

ll CRT(vector<ll> rem, vector<ll> mod, int k){
    ll m = 1;
    for(auto i : mod) m *= i;
    ll ret = 0;
    for(int i = 0 ; i < k ; i++){
        ll tmp = (m/mod[i])%mod[i];
        ll si = mod_inv(tmp,mod[i]);
        ret += (rem[i]*si%m)*(m/mod[i])%m;
        ret %= m;
    }
    return ret;
}
