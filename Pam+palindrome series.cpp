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
            insert(s[i]-'a',i);//�Ѳ����ӵ�init������ܶ�
          dp(i);
       }
    }
    int newnode(int x){
    //����һ���½ڵ㣬����Ϊx
       len[++tot]=x;return tot;
    }
    int getfail(int x,int n){
    //��failָ��֪���ҵ���׺����Ϊֹ
    while(s[n-len[x]-1]!=s[n]) x=fail[x];
      return x;
    }
    void insert(int x,int i)
    {
        //�ҵ����Ի��ĵ�λ��
        int p=getfail(last,i);
        if(!ch[p][x]){
            //�������ת�ƾͲ��ý��ˣ�����Ҫ�½�
            //ǰ�󶼼������ַ��������»��Ĵ�����Ҫ��2
            int q=newnode(len[p]+2);
            //��Ϊfailָ��ĵ���ԭ�����ϸ��׺������Ҫ��p��fail��ʼ����
            fail[q]=ch[getfail(fail[p],i)][x];
            //��¼ת��
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
            g[p][k] = f[i-len[anc[p]]-diff[p]][k]-(i-len[anc[p]]-diff[p]);//����ת���漰�µĳ���  dp[j][k-1]-j
           //    (len[anc[p]+diff[p]����Сһ�ζ����¼ӽ����Ķ�����border) Ȼ��������һ�εȲ�����ʱ ��ȥ���¼ӵ�border �ٴ���εȲ�������β��
          // ������log��
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
//����PAM��dfs
void dfs(int nl,int nr)
{
    if(nl>1&&nr>1) ans+=1ll*p1.cnt[nl]*p2.cnt[nr];//�������
    for(int i=0;i<26;i++)
    {
        if(p1.ch[nl][i]&&p2.ch[nr][i])
        dfs(p1.ch[nl][i],p2.ch[nr][i]);
    }
}
	for(int i=p1.tot;i>=2;i--) p1.cnt[p1.fail[i]]+=p1.cnt[i];//cnt��ʾ���ʲ�ͬ�Ĵ��ĸ��� 
	for(int i=p2.tot;i>=2;i--) p2.cnt[p2.fail[i]]+=p2.cnt[i];
    dfs(0,0);dfs(1,1);//������������ 
