#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int N=1e6+5;
const int mod=1e9+7;
char s[N],ss[N],s1[N];
ll ans;
struct pam
{
    int fail[N],cnt[N],len[N],tot,last,ch[N][26];
    int diff[N],anc[N],g[N][4],f[N][4];
    char s[N];
    void init(char *ss)
    {
        int n=strlen(ss+1);
        memset(fail, 0, sizeof(int) * (tot + 1));
        memset(len, 0, sizeof(int) * (tot + 1));
        memset(f, 0, sizeof(int) * (n + 1) * 4);
        memset(ch, 0, sizeof(int) * (tot + 1) * 26);
        s[0]=-1,fail[0]=1,last=0;
        len[0]=0,len[1]=-1,tot=1;
        for(int i=1;ss[i];i++)
        {
          s[i]=ss[i];        
            insert(s[i]-'a',i);//把插入扔到init里面会快很多
          dp(i);
       }
    }
    int newnode(int x){
    //建立一个新节点，长度为x
       len[++tot]=x;return tot;
    }
    int getfail(int x,int n){
    //跳fail指针知道找到后缀回文为止
    while(s[n-len[x]-1]!=s[n]) x=fail[x];
      return x;
    }
    void insert(int x,int i)
    {
        //找到可以回文的位置
        int p=getfail(last,i);
        if(!ch[p][x]){
            //如果有了转移就不用建了，否则要新建
            //前后都加上新字符，所以新回文串长度要加2
            int q=newnode(len[p]+2);
            //因为fail指向的得是原串的严格后缀，所以要从p的fail开始找起
            fail[q]=ch[getfail(fail[p],i)][x];
            //记录转移
            ch[p][x]=q;
            diff[q]=len[q]-len[fail[q]];
            anc[q]=diff[q]==diff[fail[q]]?anc[fail[q]]:fail[q];
            
        }
        cnt[last=ch[p][x]]++;
        //    printf("%d\n",last);
    }
    void dp(int i)
    {
        for(int k=1;k<=3;k++)
        {
           f[i][k]=f[i-1][k];
          for (int p=last;p>1;p =anc[p]){
            g[p][k] = f[i-len[anc[p]]-diff[p]][k]-(i-len[anc[p]]-diff[p]);//这里转移涉及新的长度  dp[j][k-1]-j
           //    (len[anc[p]+diff[p]是最小一段对于新加进来的东西的border) 然后当跳到下一段等差数列时 再去加新加的border 再从这段等差数列往尾传
          // 所以跳log次
            if (diff[p]==diff[fail[p]]){
                g[p][k]=max(g[p][k],g[fail[p]][k]);
            }
            f[i][k]=max(f[i][k],k==1?len[p]:g[p][k-1]+i);
          }    
        }
    }
}p;
int main(){
    int t;cin>>t;
    while(t--)
    {
       scanf("%s",s+1);
       p.init(s);
       int n=strlen(s+1);
       printf("%d\n",p.f[n][3]);    
     }
    return 0;
}
//两个PAM的dfs
void dfs(int nl,int nr)
{
    if(nl>1&&nr>1) ans+=1ll*p1.cnt[nl]*p2.cnt[nr];//特判虚根
    for(int i=0;i<26;i++)
    {
        if(p1.ch[nl][i]&&p2.ch[nr][i])
        dfs(p1.ch[nl][i],p2.ch[nr][i]);
    }
}
	for(int i=p1.tot;i>=2;i--) p1.cnt[p1.fail[i]]+=p1.cnt[i];//cnt表示本质不同的串的个数 
	for(int i=p2.tot;i>=2;i--) p2.cnt[p2.fail[i]]+=p2.cnt[i];
    dfs(0,0);dfs(1,1);//两个根往下找 
