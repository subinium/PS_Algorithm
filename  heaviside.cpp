// heaviside cover-up method
// (1 + a1x + a1x^2 ...)(1 + a2x + a2x^2 ...)... -> x^n coefficient?
// 1/(1-a1x) * 1/(1-a2x) ... -> q1/(1-a1x) + q2/(1-a2x) + ...
// qi = 1/(1-a2/a1) * 1/(1-a3/a1) ...
// x^n's coefficient = q1 * (a1)^n + q2 * (a2)^n ...

ll power(ll a, ll b, ll m){
    ll ret = 1;
    while(b){
        if(b&1) ret = ret*a%m;
        a = a*a%m;
        b >>= 1;
    }
    return ret;
}

ll heaviside(ll n, ll m, vector<ll> a, ll MOD){
  vector<ll> q(m);
  for(int i = 0 ; i < m ; i++){
      ll tq = 1; ll ai = power(a[i],MOD-2,MOD);
      for(int j = 0 ; j < m ; j++){
          if(i==j) continue;
          ll ttq = (1+MOD)-(ai*a[j])%MOD;
          tq *= power(ttq,MOD-2,MOD);
          tq %= MOD;
      }
      q[i] = tq;
  }
  ll ans = 0;
  for(int i = 0 ; i < m ; i++){
      ans += q[i]*power(a[i],n,MOD);
      ans %= MOD;
  }
  return ans;
}
