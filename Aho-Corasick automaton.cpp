#include<bits/stdc++.h>
using namespace std;
const int N = 300010;
const int T = 300010;
const int S = 300010;
char s[S],ss[S];
queue<int> q;
int dp[300010][20];
int n, tr[T][26], fail[T], tot, len[N], siz[T];
int ask(int l,int r)
{
    if(l>r) return 1e9+7;
    int k=log2(r-l+1);
    return min(dp[l+(1<<k)-1][k],dp[r][k]);//����
}
int main()
{
    int i,j,u;
    scanf("%d",&n);
    scanf("%s",ss+1);
    for (i=1;i<=n;i++)
    {
        scanf("%s",s);
        int l=strlen(s);
        for(u=0,j=0;s[j];++j)
        {
            int c = s[j] - 'a';
            if (!tr[u][c]) tr[u][c] = ++tot;
            u = tr[u][c];
        }
        len[u]=l;//���ַ�����β������Ϣ
    }
    for(int i=0;i<26;i++)if(tr[0][i]){fail[tr[0][i]]=0;q.push(tr[0][i]);}
    while (!q.empty())
    {
        int u = q.front();q.pop();
        len[u]=max(len[u],len[fail[u]]);
        //��Щ�㲻�ǽ�β�� ͨ����failת��len 
       //��˼���� ����A--->B���B�͸���A�ĺ��浫����Ҫת���ֲ����ַ���β ��B�� ���Ը��¹���
        for (int i=0;i<26;i++)
        {
            if (tr[u][i])
            {
                fail[tr[u][i]] = tr[fail[u]][i];
                q.push(tr[u][i]);
            }
            else tr[u][i] = tr[fail[u]][i];
        }
    }
    memset(dp,0x3f3f3f3f,sizeof(dp));dp[0][0]=0;
    int l=strlen(ss+1);u=0;
    for(int i=1;i<=l;i++)
    {
       u=tr[u][ss[i]-'a'];
       dp[i][0]=ask(i-len[u],i-1)+1;
       for(int j=1;j<20;j++)
       {
            if(i-(1<<j)+1<0) break;
            dp[i][j]=min(dp[i][j-1],dp[i-(1<<j-1)][j-1]);
        }
    }
    printf("%d",dp[l][0]>1e9?-1:dp[l][0]);
    return 0;
}
//���ַ���ST��Ҳ���Դ����߶���
