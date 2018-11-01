// tetration a^a^a ... b times mod m
const int SZ = 1e6;
const ll MOD = 1e8;

vector<ll> p;
bool ck[SZ];

void era(){
    for(int i = 2 ; i < SZ ; i++){
        if(ck[i]) continue;
        p.push_back(i);
        for(int j = 2*i ; j < SZ ; j+=i){
            ck[j] = true;
        }
    }
}

ll power(ll a, ll b, ll mod){
    ll ret = 1;
    while(b){
        if(b&1){
            ret *= a;
            if(ret > mod) ret = ret % mod + mod;
        }
        a *= a;
        if(a > mod) a = a % mod + mod;
        b >>= 1;
    }
    return ret;
}

ll toitient(ll m){
    ll ret = m;
    for(ll pr : p){
        if(pr>m) break;
        if(m%pr) continue;
        while(m%pr==0) m/=pr;
        ret = ret/pr*(pr-1);
    }
    if(m != 1) ret = ret/m*(m-1);
    return ret;
}

ll solve(ll a, ll b, ll mod){
    if(b==0) return 1;
    if(mod==1) return mod;
    ll small_ans = solve(a, b-1, toitient(mod));
    return power(a,small_ans,mod);
}
