#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
int nxt[10000005];
int s[10000005];//sΪƥ�䴮
int t[10000005];
ll n,p,ans[10000005];
int main()
{
    scanf("%lld",&n);
    for(int i=1;i<=n;i++) cin>>s[i];
    for(int i=n;i>=1;i--) t[i]=s[n-i+1];
    int k=0;//kΪ��ƥ�䵽�����ĺ�׺��ǰ׺
    nxt[1]=0;
    for(int i=2;i<=n;i++)//��ʼƥ��s
    {
        while(k!=0&&s[i]!=s[k+1])  k=nxt[k];
        if(s[i]==s[k+1])k++;
        nxt[i]=k;
    }
    k=0;
    for(int i=1;i<=n;i++)//��ʼƥ���ı���
    {
        while(k!=0&&t[i]!=s[k+1]) k=nxt[k];
        if(t[i]==s[k+1]) k++;
       if(k!=0) ans[k]++;
    }
    for(int i=n;i>=1;i--) ans[nxt[i]]+=ans[i];//�����ڵ��������
    for(int i=1;i<=n;i++) printf(i==n?"%d\n":"%d ",ans[i]);
    return 0;
}   
