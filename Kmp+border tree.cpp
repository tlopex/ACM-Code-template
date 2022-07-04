#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
int nxt[10000005];
int s[10000005];//s为匹配串
int t[10000005];
ll n,p,ans[10000005];
int main()
{
    scanf("%lld",&n);
    for(int i=1;i<=n;i++) cin>>s[i];
    for(int i=n;i>=1;i--) t[i]=s[n-i+1];
    int k=0;//k为所匹配到的最大的后缀的前缀
    nxt[1]=0;
    for(int i=2;i<=n;i++)//开始匹配s
    {
        while(k!=0&&s[i]!=s[k+1])  k=nxt[k];
        if(s[i]==s[k+1])k++;
        nxt[i]=k;
    }
    k=0;
    for(int i=1;i<=n;i++)//开始匹配文本串
    {
        while(k!=0&&t[i]!=s[k+1]) k=nxt[k];
        if(t[i]==s[k+1]) k++;
       if(k!=0) ans[k]++;
    }
    for(int i=n;i>=1;i--) ans[nxt[i]]+=ans[i];//往根节点子树求和
    for(int i=1;i<=n;i++) printf(i==n?"%d\n":"%d ",ans[i]);
    return 0;
}   
