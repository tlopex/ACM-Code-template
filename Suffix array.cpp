//���� 
#include <bits/stdc++.h>
using namespace std;
const int maxn = 2000020;
typedef long long ll;
//��������ws��������
int sa[maxn], wv[maxn], wss[maxn], wa[maxn], wb[maxn], r[maxn];
char s[maxn];
bool cmp(int *r, int a, int b, int l)
{
    return r[a] == r[b] && r[a + l] == r[b + l];
}
//O(nlogn)�����±��0��ʼ
void get_sa(int *r, int *sa, int n, int m)
{
    int *x=wa, *y=wb;
    int p =0, i, j;
    for(i = 0; i < m; i++) wss[i] = 0;
    for(i = 0; i < n; i++) wss[ x[i]=r[i] ]++;
    for(i = 1; i <= m; i++) wss[i] += wss[i - 1];
    for(i = n - 1; i >= 0; i--) sa[--wss[x[i]]] = i;
    for(j = 1, p = 1; p < n; j *= 2, m = p)
    {
        //�Եڶ��ؼ�������
        for(p = 0, i = n - j; i < n; i++) // [n-j,n)û������
            y[p++] = i;
        for(i = 0; i < n; i++)
            if(sa[i] >= j) y[p++] = sa[i] - j;
        //�Ե�һ�ؼ�������
        for(i = 0; i < n; i++) wv[i] = x[y[i]]; //����Ϊi�ĵڶ��ؼ��ֶ�Ӧ�ĵ�һ�ؼ��ֵ�������x��ʱ�൱��rk��y�൱�ڵڶ��ؼ��ֵ�sa
        for(i = 0; i < m; i++) wss[i] = 0;
        for(i = 0; i < n; i++) wss[wv[i]]++;
        for(i = 1; i <= m; i++) wss[i] += wss[i - 1];
        for(i = n - 1; i >= 0; i--) sa[--wss[wv[i]]] = y[i];
        //��ͬ���ַ���������ͬ
        swap(x,y);
        for(i = 1, p = 1, x[sa[0]] = 0; i < n; i++)
            x[sa[i]] = cmp(y, sa[i-1], sa[i], j) ? p - 1 : p++;
    }
}
//O(n)
int rk[maxn], height[maxn];
void get_height(int n)
{
    int k = 0;
    for(int i = 1; i <= n; i++) rk[sa[i]] = i;
    for(int i = 0; i < n; i++)
    {
        k ? k-- : 0;//��������height[rank[i]] �� (height[rank[i-1]] -1)
        int j = sa[rk[i] - 1];//��һ���Ŀ�ʼ�±�
        while(r[i + k] == r[j + k]) k++;
        height[rk[i]] = k;
    }
}
int main()
{
    int T;scanf("%d",&T);
    while(T--)
    {
        scanf("%s",s);
        int Max=-1;
        int len = strlen(s);
        for(int i = 0; i < len; i++)
        {
            r[i] = (int)s[i];
            if(r[i]>Max) Max=r[i];
        }
        r[len] = 0;//��Ҫ���ǣ���
        get_sa(r,sa, len +1,Max+1);// ��n+ 1
        get_height(len);
        for(int i=1;i<=len;i++) printf("%d ",sa[i]+1);
        printf("\n");
        for(int i=1;i<=len;i++) printf("%d ",height[i]);
        printf("\n");
    }
    
        
    return 0;
}

//DC3
#include<bits/stdc++.h>
using namespace std;
#define MAXLEN 4000000
#define MAXN 2000000
#define N 4000000
#define F(x) ((x)/3+((x)%3==1?0:tb))
#define G(x) ((x)<tb?(x)*3+1:((x)-tb)*3+2)
int wa[N],wb[N],wv[N],wss[N],sa[N*3];
int c0(int *r,int a,int b){return r[a]==r[b]&&r[a+1]==r[b+1]&&r[a+2]==r[b+2];}
int c12(int k,int *r,int a,int b){
    if(k==2) return r[a]<r[b]||r[a]==r[b]&&c12(1,r,a+1,b+1);
    else return r[a]<r[b]||r[a]==r[b]&&wv[a+1]<wv[b+1];
}
void sort(int *r,int *a,int *b,int n,int m)
{
    int i;
    for(int i=0;i<n;i++)wv[i]=r[a[i]];
    for(int i=0;i<m;i++)wss[i]=0;
    for(int i=0;i<n;i++)wss[wv[i]]++;
    for(int i=1;i<m;i++)wss[i]+=wss[i-1];
    for(int i=n-1;i>=0;i--)b[--wss[wv[i]]]=a[i];
}
void dc3(int *r,int *sa,int n,int m)
{
    int i,j,*rn=r+n,*san=sa+n,ta=0,tb=(n+1)/3,tbc=0,p;
    r[n]=r[n+1]=0;
    for(i=0;i<n;i++)if(i%3!=0)wa[tbc++]=i;
    sort(r+2,wa,wb,tbc,m);sort(r+1,wb,wa,tbc,m);sort(r,wa,wb,tbc,m);
    for(p=1,rn[F(wb[0])]=0,i=1;i<tbc;i++)rn[F(wb[i])]=c0(r,wb[i-1],wb[i])?p-1:p++;
    if(p<tbc)dc3(rn,san,tbc,p);else for (i=0;i<tbc;i++)san[rn[i]]=i;
    for(i=0;i<tbc;i++)if(san[i]<tb)wb[ta++]=san[i]*3;
    if(n%3==1)wb[ta++]=n-1;
    sort(r,wb,wa,ta,m);
    for(i=0;i<tbc;i++)wv[wb[i]=G(san[i])]=i;
    for(i=0,j=0,p=0;i<ta&&j<tbc;p++)sa[p]=c12(wb[j]%3,r,wa[i],wb[j])?wa[i++]:wb[j++];
    for(;i<ta;p++)sa[p]=wa[i++];
    for(;j<tbc;p++)sa[p]=wb[j++];
}

int n,len,bl[MAXLEN+10];
int s[MAXLEN+10];
int height[MAXLEN+10];
char t[MAXLEN+10];
int rk[N*3];
void Get_height(int n){
    int i,j,k=0;
    for(int i = 1; i <= n; i++) rk[sa[i]] = i;
    for(i=0;i<n;i++){
        if(!rk[i])
            height[rk[i]]=0;
        else{
            if(k)
                k--;
            for(j=sa[rk[i]-1];s[i+k]==s[j+k];k++);
            height[rk[i]]=k;
        }
    }
}

int main()
{
    int T;scanf("%d",&T);
    while(T--)
    {
      scanf("%s",t);
         len=strlen(t);
         int MAXC=-1;
         for(int i=0;i<len;i++)
        {
          s[i]=(int)t[i];
          MAXC=max(MAXC,s[i]);
        }
        s[len]=0;
        dc3(s,sa,len+1,MAXC+1);  
        Get_height(len);
        for(int i=1;i<=len;i++) printf("%d ",sa[i]+1);
        printf("\n");
        for(int i=1;i<=len;i++) printf("%d ",height[i]);
        printf("\n");    
    }
}


//RMQ
int d[maxn][20]; //2<<20����
void RMQ_ST(int n){
    for(int i=1;i<=n;i++){
        d[i][0]=height[i];
    }
    int end_j=log(n+0.0)/log(2.0);
    int end_i;
    for(int j=1;j<=end_j;j++){
        end_i=n+1-(1<<j);
        for(int i=1;i<=end_i;i++){
            d[i][j]=min(d[i][j-1],d[i+(1<<(j-1))][j-1]);
        }
    }
}
int n;
int rmq(int l,int r){//����д���е��鷳�������l��r�Ѿ��Ǻ�׺����λ��(����ȥ����rk[l],rk[r])
    if(l>r) swap(l,r);
    if(l==r) return n-sa[l]+1;//l==rһ��������������Ϊ����λ��Ϊl�ĺ�׺�����ﷵ��INFҲ�ǿ��Ե�
    l++;//������Ҫl++����Ϊ���ǵ�height[i]�Ķ���
        //������Ҫ��ѯ��λ������[l,r]֮��ĺ�׺,�����������ڵ�LCP��Сֵ��l��Ҫ++��
        //��Ϊheight[l]����λΪl�ĺ�׺����λΪl-1�ĺ�׺��LCP
        //���������ǲ���Ҫheight[l]�ģ�ֻ��Ҫheight[l+1]��height[r]
    int k=0;
    while((1<<(k+1))<=r-l+1) k++;
    int a=d[l][k];
    int b=d[r-(1<<k)+1][k];
    return min(a,b);
}
int getl(int l,int r,int len,int x)
{
    int ans=r;
    while(l<=r)
    {
        int mid=l+r>>1;
    
        if(rmq(x,mid)>=len)
        {
            ans=mid;
            r=mid-1;
        }
        else l=mid+1;
    }
    return ans;
}
int getr(int l,int r,int len,int x)
{
    int ans=l;

    while(l<=r)
    {
        int mid=(l+r)>>1;
    
        if(rmq(x,mid)>=len)
        {
            ans=mid;
            l=mid+1;
        }
        else r=mid-1;
    }
    return ans;
}
