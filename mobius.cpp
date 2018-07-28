#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
typedef long long ll;

int mu[1000001];
ll a[1000001];

int main(){
    mu[1] = 1;
    for(int i = 1;i<=1000000;i++){
        for(int j = 2*i;j<=1000000;j += i){
            mu[j] -= mu[i];
        }
    }
    int n;
    scanf("%d",&n);
    for(int i = 0,tmp ; i < n ; i++){
        scanf("%d",&tmp);
        for(int j = 1 ; j*j <= tmp ; j++){
            if(tmp%j==0){
                a[j]++;
                a[tmp/j]++;
            }
            if(tmp==j*j) a[j]--;
        }
    }
    ll tot = 0;
    for(int i = 1 ; i <= 1000000 ; i++){
        if(a[i]>2) tot += (ll)mu[i]*(a[i])*(a[i]-1)*(a[i]-2)/6ll;
    }
    printf("%lld",tot);
}
