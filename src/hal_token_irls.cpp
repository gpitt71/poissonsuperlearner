// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#include <cmath>
using namespace Rcpp;

namespace {
struct Fenwick {
  std::vector<double> a,b;
  explicit Fenwick(int n):a(n+1,0),b(n+1,0){}
  void add(int ix,double x,double y){for(int k=ix;k<(int)a.size();k+=k&-k){a[k]+=x;b[k]+=y;}}
  std::pair<double,double> sum(int ix) const {double x=0,y=0;for(int k=ix;k>0;k-=k&-k){x+=a[k];y+=b[k];}return {x,y};}
};
struct Token{int var;int kind;double cut;int level;};
bool passes(const Token& t, const NumericMatrix& x, const IntegerMatrix& cat, int i) {
  return t.kind==0 ? (std::isfinite(x(i,t.var)) && x(i,t.var)>=t.cut) :
                     (cat(i,t.var)!=NA_INTEGER && cat(i,t.var)==t.level);
}
std::vector<int> eligible(const Token& t,const NumericMatrix& x,const IntegerMatrix& cat,int n) {
  std::vector<int> ids;ids.reserve(n);
  for(int i=0;i<n;++i) if(passes(t,x,cat,i)) ids.push_back(i);
  return ids;
}
}

// Numeric variables reside in numeric_data; categorical variables in cat_data.
// var ids are 0-based shared column indexes; unused entries in either matrix are ignored.
// A basis column is one or two tokens, supplied by token metadata vectors.
// [[Rcpp::export]]
List hal_token_sums_cpp(NumericMatrix numeric_data, IntegerMatrix cat_data,
    IntegerVector var1, IntegerVector kind1, NumericVector cut1, IntegerVector level1,
    IntegerVector var2, IntegerVector kind2, NumericVector cut2, IntegerVector level2,
    NumericVector weights, NumericVector signal) {
  int n=numeric_data.nrow(), d=numeric_data.ncol(), p=var1.size();
  if(cat_data.nrow()!=n||cat_data.ncol()!=d||weights.size()!=n||signal.size()!=n)
    stop("Mismatched data dimensions.");
  if(kind1.size()!=p||cut1.size()!=p||level1.size()!=p||var2.size()!=p||kind2.size()!=p||cut2.size()!=p||level2.size()!=p)
    stop("Mismatched token lengths.");
  std::vector<Token> first(p),second(p);
  for(int j=0;j<p;++j){
    if(var1[j]<0||var1[j]>=d||kind1[j]<0||kind1[j]>1)stop("Invalid first token.");
    first[j]={var1[j],kind1[j],cut1[j],level1[j]};
    if(var2[j]>=0){
      if(var2[j]>=d||var2[j]==var1[j]||kind2[j]<0||kind2[j]>1)stop("Invalid second token.");
      second[j]={var2[j],kind2[j],cut2[j],level2[j]};
    } else second[j]={-1,-1,0,-1};
  }
  for(int i=0;i<n;++i)if(!std::isfinite(weights[i])||!std::isfinite(signal[i]))stop("Nonfinite input weights or signal.");
  NumericVector curvature(p),gradient(p);
  // Group columns by their ordered variable/kind tuple. No HAL columns constructed.
  std::unordered_map<std::string,std::vector<int>> groups;
  for(int j=0;j<p;++j){
    std::string key=std::to_string(first[j].var)+":"+std::to_string(first[j].kind)+"/"+
      std::to_string(second[j].var)+":"+std::to_string(second[j].kind);
    groups[key].push_back(j);
  }
  std::vector<int> rows(n); std::iota(rows.begin(),rows.end(),0);
  for(auto &entry:groups){
    const auto &js=entry.second;
    const Token a=first[js[0]],b=second[js[0]];
    if(b.var<0){
      if(a.kind==0){
        // Sort once per group and sweep thresholds from largest to smallest.
        std::sort(rows.begin(),rows.end(),[&](int i,int k){return numeric_data(i,a.var)>numeric_data(k,a.var);});
        std::vector<int> order=js;
        std::sort(order.begin(),order.end(),[&](int i,int k){return first[i].cut>first[k].cut;});
        int at=0;double wsum=0,ssum=0;
        for(int j:order){while(at<n && numeric_data(rows[at],a.var)>=first[j].cut){
          int i=rows[at++];wsum+=weights[i];ssum+=signal[i];
        }curvature[j]=wsum/n;gradient[j]=ssum/n;}
      }else{
        std::unordered_map<int,std::pair<double,double>> sums;
        for(int i=0;i<n;++i){int lev=cat_data(i,a.var);if(lev==NA_INTEGER)continue;
          sums[lev].first+=weights[i];sums[lev].second+=signal[i];}
        for(int j:js){auto it=sums.find(first[j].level);if(it!=sums.end()){
          curvature[j]=it->second.first/n;gradient[j]=it->second.second/n;
        }}
      }
      continue;
    }
    if(a.kind==1&&b.kind==1){
      // Category x category: aggregate each observed pair once.
      std::unordered_map<long long,std::pair<double,double>> sums;
      for(int i=0;i<n;++i){int av=cat_data(i,a.var),bv=cat_data(i,b.var);
        if(av==NA_INTEGER||bv==NA_INTEGER)continue;
        long long key=((long long)av<<32)|(unsigned int)bv;
        sums[key].first+=weights[i];sums[key].second+=signal[i];}
      for(int j:js){long long key=((long long)first[j].level<<32)|(unsigned int)second[j].level;
        auto it=sums.find(key);if(it!=sums.end()){
          curvature[j]=it->second.first/n;gradient[j]=it->second.second/n;
        }}
      continue;
    }
    if(a.kind!=b.kind){
      // Partition by categorical level; inside each partition traverse numeric
      // observations and cutpoints from largest to smallest.
      bool first_numeric=a.kind==0;
      int nc=first_numeric?a.var:b.var, cc=first_numeric?b.var:a.var;
      std::unordered_map<int,std::vector<int>> bycat, bycol;
      for(int i=0;i<n;++i){int lev=cat_data(i,cc);if(lev!=NA_INTEGER)bycat[lev].push_back(i);}
      for(int j:js)bycol[first_numeric?second[j].level:first[j].level].push_back(j);
      for(auto &part:bycol){
        auto it=bycat.find(part.first);if(it==bycat.end())continue;
        auto &r=it->second;
        std::sort(r.begin(),r.end(),[&](int i,int k){return numeric_data(i,nc)>numeric_data(k,nc);});
        auto &cols=part.second;
        std::sort(cols.begin(),cols.end(),[&](int i,int k){return (first_numeric?first[i].cut:second[i].cut)>(first_numeric?first[k].cut:second[k].cut);});
        size_t at=0;double wsum=0,ssum=0;
        for(int j:cols){double cut=first_numeric?first[j].cut:second[j].cut;
          while(at<r.size()&&numeric_data(r[at],nc)>=cut){int i=r[at++];wsum+=weights[i];ssum+=signal[i];}
          curvature[j]=wsum/n;gradient[j]=ssum/n;
        }
      }
      continue;
    }
    // Continuous x continuous: sweep descending X1 and add new rows to a
    // Fenwick tree indexed by ascending X2 values; query suffix mass for X2>=cut.
    std::vector<int> finite;finite.reserve(n);
    for(int i=0;i<n;++i)if(std::isfinite(numeric_data(i,a.var))&&std::isfinite(numeric_data(i,b.var)))finite.push_back(i);
    std::vector<double> values2;values2.reserve(finite.size());
    for(int i:finite)values2.push_back(numeric_data(i,b.var));
    std::sort(values2.begin(),values2.end());values2.erase(std::unique(values2.begin(),values2.end()),values2.end());
    std::sort(finite.begin(),finite.end(),[&](int i,int k){return numeric_data(i,a.var)>numeric_data(k,a.var);});
    std::vector<int> cols=js;
    std::sort(cols.begin(),cols.end(),[&](int i,int k){return first[i].cut>first[k].cut;});
    Fenwick bit(values2.size());size_t at=0;double totalw=0,totals=0;
    for(int j:cols){
      while(at<finite.size()&&numeric_data(finite[at],a.var)>=first[j].cut){
        int i=finite[at++];int ix=std::lower_bound(values2.begin(),values2.end(),numeric_data(i,b.var))-values2.begin()+1;
        bit.add(ix,weights[i],signal[i]);totalw+=weights[i];totals+=signal[i];
      }
      int before=std::lower_bound(values2.begin(),values2.end(),second[j].cut)-values2.begin();
      auto lower=bit.sum(before);
      curvature[j]=(totalw-lower.first)/n;gradient[j]=(totals-lower.second)/n;
    }
  }
  return List::create(_["curvature"]=curvature,_["gradient"]=gradient);
}

// A direct token support oracle for validating incremental residual updates;
// deliberately O(n * interaction_degree) and not a fast coordinate-descent kernel.
// [[Rcpp::export]]
IntegerVector hal_token_support_cpp(NumericMatrix numeric_data,IntegerMatrix cat_data,
    int var1,int kind1,double cut1,int level1,int var2=-1,int kind2=-1,double cut2=0,int level2=-1){
  int n=numeric_data.nrow();Token a{var1,kind1,cut1,level1},b{var2,kind2,cut2,level2};
  if(cat_data.nrow()!=n||cat_data.ncol()!=numeric_data.ncol())stop("Dimensions mismatch");
  if(var1<0||var1>=numeric_data.ncol()||(var2>=0&&(var2>=numeric_data.ncol()||var2==var1)))stop("Invalid variable ids");
  std::vector<int> out;
  for(int i=0;i<n;++i)if(passes(a,numeric_data,cat_data,i)&&(var2<0||passes(b,numeric_data,cat_data,i)))out.push_back(i+1);
  return wrap(out);
}

// Cached ordered indices: construct once, reuse at every IRLS evaluation.
#include <memory>
#include <map>
#include <utility>
namespace {
struct PairIndex {
 int va,vb;
 std::vector<int> rows, ranks;
 std::vector<double> values;
};
struct TripleIndex {
 int vars[3]{-1,-1,-1}, kinds[3]{-1,-1,-1}, dims[3]{0,0,0};
 std::vector<double> cuts[3];
 std::vector<int> levels[3];
 // Static row/column positions in a degree-three threshold grid. Dynamic
 // weights are accumulated into this grid by each derivative evaluation.
 std::vector<int> row_cell, col_cell;
};
struct Prepared {
 int n,p;
 std::vector<std::vector<Token>> tokens;
 std::vector<std::vector<int>> num_sorted;
 std::vector<std::map<int,std::vector<int>>> categories;
 std::map<std::pair<int,int>,PairIndex> pairs;
 std::vector<std::vector<int>> group_cols;
 std::vector<int> group_type, group_v1, group_v2;
 std::vector<TripleIndex> triples;
 std::vector<std::map<int,std::vector<int>>> mixed_rows;
 std::vector<std::map<int,std::vector<int>>> mixed_cols;
};
}

// [[Rcpp::export]]
SEXP hal_token_index_cpp(NumericMatrix numeric_data, IntegerMatrix cat_data,
 IntegerVector token_offsets,IntegerVector token_var,IntegerVector token_kind,
 NumericVector token_cut,IntegerVector token_level){
 const int n=numeric_data.nrow(),d=numeric_data.ncol(),p=token_offsets.size()-1;
 const int nt=token_var.size();
 if(cat_data.nrow()!=n||cat_data.ncol()!=d||p<1||token_offsets[0]!=0||
    token_offsets[p]!=nt||token_kind.size()!=nt||token_cut.size()!=nt||
    token_level.size()!=nt)stop("Invalid token dimensions");
 auto *P=new Prepared();
 XPtr<Prepared> ptr(P,true);
 P->n=n;P->p=p;P->tokens.resize(p);P->num_sorted.resize(d);P->categories.resize(d);
 std::map<std::string,int> keys;
 for(int j=0;j<p;++j){
  const int begin=token_offsets[j],end=token_offsets[j+1];
  if(end<=begin)stop("Empty token column");
  std::vector<unsigned char> seen(d,0);
  for(int k=begin;k<end;++k){
   if(token_var[k]<0||token_var[k]>=d||token_kind[k]<0||token_kind[k]>1||seen[token_var[k]])
    stop("Invalid token");
   seen[token_var[k]]=1;
   P->tokens[j].push_back({token_var[k],token_kind[k],token_cut[k],token_level[k]});
  }
  const Token &a=P->tokens[j][0];
  const bool has_second=P->tokens[j].size()>1;
  const Token b=has_second?P->tokens[j][1]:Token{-1,-1,0,-1};
  // Keep degree-three variable/kind combinations separate so that their
  // complete threshold grids can be evaluated collectively.
  std::string key=std::to_string(P->tokens[j].size());
  if(P->tokens[j].size()<=3){
   for(const Token &token:P->tokens[j])
    key += ":"+std::to_string(token.var)+":"+std::to_string(token.kind);
  }else key += ":fallback";
  auto it=keys.find(key);
  if(it==keys.end()){
   int g=P->group_cols.size();keys[key]=g;P->group_cols.emplace_back();
   P->group_type.push_back(P->tokens[j].size()==3?5:(P->tokens[j].size()>3?6:
    (!has_second?(a.kind==0?0:1):(a.kind+b.kind==0?4:(a.kind+b.kind==2?2:3)))));
   P->group_v1.push_back(a.var);P->group_v2.push_back(b.var);
   P->mixed_rows.emplace_back();P->mixed_cols.emplace_back();it=keys.find(key);
  }
  P->group_cols[it->second].push_back(j);
 }
 for(int v=0;v<d;++v){
  auto &ord=P->num_sorted[v];ord.reserve(n);
  for(int i=0;i<n;++i){
   double x=numeric_data(i,v);
   if(std::isfinite(x))ord.push_back(i);
   const int c=cat_data(i,v);if(c!=NA_INTEGER)P->categories[v][c].push_back(i);
  }
  std::sort(ord.begin(),ord.end(),[&](int i,int k){return numeric_data(i,v)>numeric_data(k,v);});
 }
 // Cache numeric descending rows per categorical level for mixed interaction groups.
 for(size_t g=0;g<P->group_cols.size();++g)if(P->group_type[g]==3){
  int j=P->group_cols[g][0];const Token &a=P->tokens[j][0],&b=P->tokens[j][1];
  int nv=a.kind==0?a.var:b.var;int cv=a.kind==1?a.var:b.var;
  for(int i:P->num_sorted[nv]){int c=cat_data(i,cv);if(c!=NA_INTEGER)P->mixed_rows[g][c].push_back(i);}
  const bool first_num=a.kind==0;
  for(int col:P->group_cols[g])P->mixed_cols[g][first_num?P->tokens[col][1].level:P->tokens[col][0].level].push_back(col);
  for(auto &entry:P->mixed_cols[g])std::sort(entry.second.begin(),entry.second.end(),[&](int x,int y){
   return (first_num?P->tokens[x][0].cut:P->tokens[x][1].cut)>
          (first_num?P->tokens[y][0].cut:P->tokens[y][1].cut);});
 }
 // Store one sorted row order and X2 coordinate ranks per continuous pair.
 for(size_t g=0;g<P->group_cols.size();++g) if(P->group_type[g]==4){
  int va=P->group_v1[g],vb=P->group_v2[g];auto key=std::make_pair(va,vb);
  if(P->pairs.count(key))continue;
  PairIndex pi;pi.va=va;pi.vb=vb;
  for(int i:P->num_sorted[va])if(std::isfinite(numeric_data(i,vb))){
   pi.rows.push_back(i);pi.values.push_back(numeric_data(i,vb));
  }
  std::sort(pi.values.begin(),pi.values.end());pi.values.erase(std::unique(pi.values.begin(),pi.values.end()),pi.values.end());
  pi.ranks.reserve(pi.rows.size());
  for(int i:pi.rows)pi.ranks.push_back(1+std::lower_bound(pi.values.begin(),pi.values.end(),numeric_data(i,vb))-pi.values.begin());
  P->pairs.emplace(key,std::move(pi));
 }
 // Threshold order never changes across IRLS iterations.
 for(size_t g=0;g<P->group_cols.size();++g)if(P->group_type[g]==0||P->group_type[g]==4){
  auto &cols=P->group_cols[g];std::sort(cols.begin(),cols.end(),[&](int x,int y){
   return P->tokens[x][0].cut>P->tokens[y][0].cut;});
 }
 // Degree-three columns form a complete grid within each variable/kind
 // combination. Cache grid coordinates for rows and columns once. Numeric
 // axes will be suffix-summed; categorical axes remain exact partitions.
 P->triples.resize(P->group_cols.size());
 for(size_t g=0;g<P->group_cols.size();++g)if(P->group_type[g]==5){
  TripleIndex &ti=P->triples[g];const auto &js=P->group_cols[g];
  for(int axis=0;axis<3;++axis){
   const Token &token=P->tokens[js[0]][axis];
   ti.vars[axis]=token.var;ti.kinds[axis]=token.kind;
   if(token.kind==0){
    for(int j:js)ti.cuts[axis].push_back(P->tokens[j][axis].cut);
    auto &z=ti.cuts[axis];std::sort(z.begin(),z.end());
    z.erase(std::unique(z.begin(),z.end()),z.end());ti.dims[axis]=z.size();
   }else{
    for(int j:js)ti.levels[axis].push_back(P->tokens[j][axis].level);
    auto &z=ti.levels[axis];std::sort(z.begin(),z.end());
    z.erase(std::unique(z.begin(),z.end()),z.end());ti.dims[axis]=z.size();
   }
  }
  const size_t cells=static_cast<size_t>(ti.dims[0])*ti.dims[1]*ti.dims[2];
  if(cells>static_cast<size_t>(std::numeric_limits<int>::max()))
   stop("Degree-three token grid is too large");
  ti.col_cell.reserve(js.size());
  for(int j:js){int coord[3];
   for(int axis=0;axis<3;++axis){const Token &token=P->tokens[j][axis];
    coord[axis]=token.kind==0 ?
     std::lower_bound(ti.cuts[axis].begin(),ti.cuts[axis].end(),token.cut)-ti.cuts[axis].begin() :
     std::lower_bound(ti.levels[axis].begin(),ti.levels[axis].end(),token.level)-ti.levels[axis].begin();
   }
   ti.col_cell.push_back((coord[0]*ti.dims[1]+coord[1])*ti.dims[2]+coord[2]);
  }
  ti.row_cell.assign(n,-1);
  for(int i=0;i<n;++i){int coord[3];bool ok=true;
   for(int axis=0;axis<3&&ok;++axis){
    if(ti.kinds[axis]==0){double value=numeric_data(i,ti.vars[axis]);
     if(!std::isfinite(value)){ok=false;break;}
     coord[axis]=std::upper_bound(ti.cuts[axis].begin(),ti.cuts[axis].end(),value)-ti.cuts[axis].begin()-1;
     if(coord[axis]<0)ok=false;
    }else{int value=cat_data(i,ti.vars[axis]);
     auto at=std::lower_bound(ti.levels[axis].begin(),ti.levels[axis].end(),value);
     if(value==NA_INTEGER||at==ti.levels[axis].end()||*at!=value)ok=false;
     else coord[axis]=at-ti.levels[axis].begin();
    }
   }
   if(ok)ti.row_cell[i]=(coord[0]*ti.dims[1]+coord[1])*ti.dims[2]+coord[2];
  }
 }
 return ptr;
}

// Native bytes retained by the external pointer (vector capacities included).
// [[Rcpp::export]]
double hal_token_index_bytes_cpp(SEXP index){
 XPtr<Prepared> P(index);size_t bytes=sizeof(Prepared);
 for(const auto &z:P->tokens)bytes+=sizeof(std::vector<Token>)+z.capacity()*sizeof(Token);
 for(const auto &z:P->num_sorted)bytes+=sizeof(std::vector<int>)+z.capacity()*sizeof(int);
 for(const auto &m:P->categories)for(const auto &z:m)bytes+=sizeof(z)+z.second.capacity()*sizeof(int);
 for(const auto &z:P->pairs)bytes+=sizeof(z)+z.second.rows.capacity()*sizeof(int)+
   z.second.ranks.capacity()*sizeof(int)+z.second.values.capacity()*sizeof(double);
 for(const auto &z:P->triples){
  bytes+=sizeof(TripleIndex)+z.row_cell.capacity()*sizeof(int)+z.col_cell.capacity()*sizeof(int);
  for(int axis=0;axis<3;++axis)
   bytes+=z.cuts[axis].capacity()*sizeof(double)+z.levels[axis].capacity()*sizeof(int);
 }
 for(const auto &z:P->group_cols)bytes+=sizeof(std::vector<int>)+z.capacity()*sizeof(int);
 for(const auto &m:P->mixed_rows)for(const auto &z:m)bytes+=sizeof(z)+z.second.capacity()*sizeof(int);
 for(const auto &m:P->mixed_cols)for(const auto &z:m)bytes+=sizeof(z)+z.second.capacity()*sizeof(int);
 return static_cast<double>(bytes);
}

// [[Rcpp::export]]
List hal_token_sums_cached_cpp(SEXP index,NumericMatrix numeric_data,IntegerMatrix cat_data,
 NumericVector weights,NumericVector signal){
 XPtr<Prepared> P(index);const int n=P->n,p=P->p;
 if(weights.size()!=n||signal.size()!=n||numeric_data.nrow()!=n||cat_data.nrow()!=n)stop("Length mismatch");
 for(int i=0;i<n;++i)if(!std::isfinite(weights[i])||!std::isfinite(signal[i]))stop("Nonfinite weights/signal");
 NumericVector curvature(p),gradient(p);const double inv=1.0/n;
 for(size_t g=0;g<P->group_cols.size();++g){
  const auto &js=P->group_cols[g];int type=P->group_type[g],va=P->group_v1[g],vb=P->group_v2[g];
  if(type==0){
   const auto &rows=P->num_sorted[va];size_t at=0;double w=0,s=0;
   for(int j:js){while(at<rows.size()&&numeric_data(rows[at],va)>=P->tokens[j][0].cut){int i=rows[at++];w+=weights[i];s+=signal[i];}curvature[j]=w*inv;gradient[j]=s*inv;}
  }else if(type==1){
   for(int j:js){double w=0,s=0;auto it=P->categories[va].find(P->tokens[j][0].level);
    if(it!=P->categories[va].end())for(int i:it->second){w+=weights[i];s+=signal[i];}
    curvature[j]=w*inv;gradient[j]=s*inv;}
  }else if(type==2){
   std::map<std::pair<int,int>,std::pair<double,double>> sums;
   for(int i=0;i<n;++i){int a=cat_data(i,va),b=cat_data(i,vb);if(a!=NA_INTEGER&&b!=NA_INTEGER){
    auto &z=sums[std::make_pair(a,b)];z.first+=weights[i];z.second+=signal[i];}}
   for(int j:js){auto it=sums.find(std::make_pair(P->tokens[j][0].level,P->tokens[j][1].level));
    if(it!=sums.end()){curvature[j]=it->second.first*inv;gradient[j]=it->second.second*inv;}}
  }else if(type==3){
   const bool first_num=P->tokens[js[0]][0].kind==0;int nv=first_num?va:vb;
   // Iterate the pre-sorted numeric ordering, selecting only the category of interest.
   for(const auto &entry:P->mixed_cols[g]){int category=entry.first;const auto &cols=entry.second;
    auto found=P->mixed_rows[g].find(category);
    if(found==P->mixed_rows[g].end())continue;
    const auto &rows=found->second;
    size_t at=0;double w=0,s=0;
    for(int j:cols){double cut=first_num?P->tokens[j][0].cut:P->tokens[j][1].cut;
     while(at<rows.size()&&numeric_data(rows[at],nv)>=cut){int i=rows[at++];w+=weights[i];s+=signal[i];}
     curvature[j]=w*inv;gradient[j]=s*inv;}
   }
  }else if(type==4){
   const auto &pi=P->pairs.at(std::make_pair(va,vb));Fenwick bit(pi.values.size());
   size_t at=0;double w=0,s=0;
   for(int j:js){while(at<pi.rows.size()&&numeric_data(pi.rows[at],va)>=P->tokens[j][0].cut){
     int i=pi.rows[at];bit.add(pi.ranks[at],weights[i],signal[i]);w+=weights[i];s+=signal[i];++at;
    }
    int before=std::lower_bound(pi.values.begin(),pi.values.end(),P->tokens[j][1].cut)-pi.values.begin();
    auto less=bit.sum(before);curvature[j]=(w-less.first)*inv;gradient[j]=(s-less.second)*inv;
   }
  }else if(type==5){
   const TripleIndex &ti=P->triples[g];
   const int d0=ti.dims[0],d1=ti.dims[1],d2=ti.dims[2];
   const size_t cells=static_cast<size_t>(d0)*d1*d2;
   std::vector<double> wsum(cells,0.0),ssum(cells,0.0);
   for(int i=0;i<n;++i){int cell=ti.row_cell[i];if(cell>=0){
    wsum[cell]+=weights[i];ssum[cell]+=signal[i];
   }}
   if(ti.kinds[0]==0)for(int a=d0-2;a>=0;--a)for(int b=0;b<d1;++b)for(int c=0;c<d2;++c){
    size_t here=(static_cast<size_t>(a)*d1+b)*d2+c;
    size_t next=here+static_cast<size_t>(d1)*d2;
    wsum[here]+=wsum[next];ssum[here]+=ssum[next];
   }
   if(ti.kinds[1]==0)for(int a=0;a<d0;++a)for(int b=d1-2;b>=0;--b)for(int c=0;c<d2;++c){
    size_t here=(static_cast<size_t>(a)*d1+b)*d2+c,next=here+d2;
    wsum[here]+=wsum[next];ssum[here]+=ssum[next];
   }
   if(ti.kinds[2]==0)for(int a=0;a<d0;++a)for(int b=0;b<d1;++b)for(int c=d2-2;c>=0;--c){
    size_t here=(static_cast<size_t>(a)*d1+b)*d2+c,next=here+1;
    wsum[here]+=wsum[next];ssum[here]+=ssum[next];
   }
   for(size_t k=0;k<js.size();++k){int j=js[k],cell=ti.col_cell[k];
    curvature[j]=wsum[cell]*inv;gradient[j]=ssum[cell]*inv;
   }
  }else{
   // Exact fallback for interaction degrees above three.
   for(int j:js){double w=0,s=0;for(int i=0;i<n;++i){bool ok=true;
    for(const Token &token:P->tokens[j])if(!passes(token,numeric_data,cat_data,i)){ok=false;break;}
    if(ok){w+=weights[i];s+=signal[i];}}
    curvature[j]=w*inv;gradient[j]=s*inv;}
  }
 }
 return List::create(_["curvature"]=curvature,_["gradient"]=gradient);
}

// Token-only Poisson IRLS.  The cached Prepared index is retained for the path.
#include <chrono>
namespace {
inline double st(double x,double t){return x>t?x-t:(x< -t?x+t:0.0);}
inline std::vector<int> token_rows(const Prepared& P,const NumericMatrix& x,const IntegerMatrix& cat,int j){
 const auto &tokens=P.tokens[j];const std::vector<int>* base=nullptr;
 for(const Token &token:tokens){
  const std::vector<int>* candidate=nullptr;
  if(token.kind==0)candidate=&P.num_sorted[token.var];
  else {auto it=P.categories[token.var].find(token.level);
   if(it==P.categories[token.var].end())return {};candidate=&it->second;}
  if(base==nullptr||candidate->size()<base->size())base=candidate;
 }
 std::vector<int> out;out.reserve(base->size());
 for(int i:*base){bool ok=true;for(const Token &token:tokens)if(!passes(token,x,cat,i)){ok=false;break;}
  if(ok)out.push_back(i);}
 return out;
}
}

// [[Rcpp::export]]
List hal_token_irls_cpp(SEXP index,NumericMatrix numeric_data,IntegerMatrix cat_data,
  NumericVector y,NumericVector offset,NumericVector lambda,NumericVector penalty_factor,
  bool intercept=true,double tol=1e-8,int maxit_irls=100,int maxit_cd=10000,
  bool active_set=true,bool rescale_penalty_factor=true,
  bool strong_rules=true){
 XPtr<Prepared> P(index);const int n=P->n,p=P->p;
 if(n<1||p<1||numeric_data.nrow()!=n||cat_data.nrow()!=n||y.size()!=n||offset.size()!=n||penalty_factor.size()!=p||lambda.size()<1)stop("Invalid dimensions");
 if(!(tol>0)||maxit_irls<1||maxit_cd<1)stop("Invalid iteration controls");
 double pfs=0,ys=0,exposure=0;
 for(int i=0;i<n;++i){if(!std::isfinite(y[i])||y[i]<0||!std::isfinite(offset[i])||offset[i]>700)stop("Invalid response or offset");ys+=y[i];exposure+=std::exp(offset[i]);}
 if(intercept&&ys<=0)stop("No finite intercept for all-zero response");
 for(int j=0;j<p;++j){if(!std::isfinite(penalty_factor[j])||penalty_factor[j]<0)stop("Invalid penalty factor");pfs+=penalty_factor[j];}
 if(rescale_penalty_factor&&pfs<=0)stop("All penalties are zero");
 const double scale=rescale_penalty_factor?p/pfs:1.0;
 std::vector<double> pf(p);for(int j=0;j<p;++j)pf[j]=scale*penalty_factor[j];
 for(int k=0;k<lambda.size();++k)if(!std::isfinite(lambda[k])||lambda[k]<0||(k>0&&lambda[k]>lambda[k-1]))stop("Lambda must be finite, nonnegative and decreasing");
 NumericMatrix B(p,lambda.size());NumericVector a0(lambda.size()),objective_path(lambda.size()),loss_path(lambda.size()),pen_path(lambda.size()),kkt_path(lambda.size());
 IntegerVector iterations(lambda.size()),sweeps(lambda.size()),active_size(lambda.size());LogicalVector converged(lambda.size()),cd_converged(lambda.size());
 NumericMatrix timing(lambda.size(),10);colnames(timing)=CharacterVector::create("irls_update","curvature","coordinate_gradients","residual_updates","active_and_convergence_checks","total","support_construction","inactive_kkt","final_kkt","likelihood");
 NumericMatrix counts(lambda.size(),9);colnames(counts)=CharacterVector::create("support_builds","support_cache_hits","support_rows_scanned","coordinate_visits","coordinate_rows_scanned","inactive_kkt_scans","active_support_entries_peak","cache_rejections","screening_scans");
 IntegerVector strong_candidates(lambda.size()),strong_excluded(lambda.size()),strong_kkt_reentries(lambda.size()),strong_rule_scans(lambda.size());
 double support_ms=0; long long support_builds=0,support_hits=0,support_rows=0;
 const double inv=1.0/n;
 std::vector<double> beta(p,0),eta(n),mu(n),r(n),wr(n),zeros(n,0),ones(n,1);
 double b0=intercept?std::log(ys/exposure):0;
 std::vector<std::vector<int>> supports(p);
 std::vector<unsigned char> support_ready(p,0);
 std::vector<int> temporary_support;
 // Exact deterministic active support cache; no user-tuned capacity or eviction.
 long long cached_entries=0,peak_cached_entries=0;
 auto support=[&](int j)->const std::vector<int>&{
   if(support_ready[j]) {++support_hits; return supports[j];}
   auto t0=std::chrono::steady_clock::now();
   temporary_support=token_rows(*P,numeric_data,cat_data,j);
   support_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t0).count();
   ++support_builds; support_rows+=temporary_support.size();
   return temporary_support;
 };
 auto promote=[&](int j,const std::vector<int>& rows){
   if(support_ready[j]) return;
   supports[j]=rows;support_ready[j]=1;
   cached_entries+=rows.size();
   peak_cached_entries=std::max(peak_cached_entries,cached_entries);
  };
 auto full_sums=[&](const std::vector<double>& w,const std::vector<double>& s)->List{
  return hal_token_sums_cached_cpp(index,numeric_data,cat_data,wrap(w),wrap(s));};
 auto likelihood=[&](double lam,bool fill)->double{
  std::fill(eta.begin(),eta.end(),b0);
  for(int i=0;i<n;++i)eta[i]+=offset[i];
  for(int j=0;j<p;++j)if(beta[j]!=0){for(int i:support(j))eta[i]+=beta[j];}
  double val=0;for(int i=0;i<n;++i){if(!std::isfinite(eta[i])||eta[i]>700)return R_PosInf;
    double m=std::exp(std::max(-745.0,eta[i]));if(fill)mu[i]=m;val+=m-y[i]*eta[i];}
  for(int j=0;j<p;++j)val+=n*lam*pf[j]*std::abs(beta[j]);return val*inv;
 };
 auto kkt=[&](double lam)->double{
  std::vector<double> signal(n);for(int i=0;i<n;++i)signal[i]=mu[i]-y[i];
  auto sums=full_sums(ones,signal);NumericVector gradient=sums["gradient"];
  double mx=0;for(int j=0;j<p;++j){double g=gradient[j],bound=lam*pf[j],v;
   if(beta[j]>1e-12)v=std::abs(g+bound);else if(beta[j]<-1e-12)v=std::abs(g-bound);else v=std::max(0.0,std::abs(g)-bound);
   mx=std::max(mx,v);}
  if(intercept){double g=0;for(int i=0;i<n;++i)g+=signal[i];mx=std::max(mx,std::abs(g*inv));}return mx;
 };
 using Clock2=std::chrono::steady_clock;
 auto milliseconds=[](Clock2::time_point start){return std::chrono::duration<double,std::milli>(Clock2::now()-start).count();};
 std::vector<double> preceding_gradient(p,0.0);
 for(int kk=0;kk<lambda.size();++kk){
  auto full_start=Clock2::now();double tm_irls=0,tm_curv=0,tm_grad=0,tm_res=0,tm_checks=0,tm_inactive=0,tm_final_kkt=0,tm_likelihood=0;
  auto support_ms_begin=support_ms;long long builds_begin=support_builds,hits_begin=support_hits,rows_begin=support_rows;long long visits=0,coordinate_rows=0,inactive_scans=0,screening_scans=0;
  const double lam=lambda[kk];bool outer_ok=false,all_inner=true;int total_sweeps=0;double final_kkt=R_PosInf;
  // Sequential strong rule from the previous, fully KKT-checked solution.
  //  |g_j(lambda[k-1])| < pf[j] * (2*lambda[k] - lambda[k-1])
  // is a heuristic discard, NOT a safe rule: a full Poisson KKT check later
  // reactivates every erroneously discarded coordinate.
  std::vector<unsigned char> eligible(p,1);
  if(strong_rules && active_set && kk>0 && converged[kk-1]){
   double cutoff=2*lam-lambda[kk-1];
   if(cutoff>0){
    for(int j=0;j<p;++j){
     if(pf[j]>0 && beta[j]==0 && std::abs(preceding_gradient[j])<pf[j]*cutoff)
       eligible[j]=0;
    }
   }
  }
  auto count_eligible=[&](){int total=0;for(unsigned char v:eligible)total+=v;return total;};
  strong_candidates[kk]=count_eligible();strong_excluded[kk]=p-strong_candidates[kk];
  int reentries=0;

  for(int outer=0;outer<maxit_irls;++outer){
   if(outer%10==0)checkUserInterrupt();
   std::vector<double> old=beta;double old0=b0;
   auto tic=Clock2::now();double current=likelihood(lam,true);tm_likelihood+=milliseconds(tic);
   if(!std::isfinite(current))stop("Nonfinite Poisson objective");
   for(int i=0;i<n;++i){if(mu[i]==0&&y[i]>0)stop("Underflowed positive-response mean");r[i]=mu[i]>0?(y[i]-mu[i])/mu[i]:-1;}
   tm_irls+=milliseconds(tic);
   tic=Clock2::now();NumericVector curvature=full_sums(mu,zeros)["curvature"];tm_curv+=milliseconds(tic);
   std::vector<unsigned char> active(p);for(int j=0;j<p;++j)active[j]=(beta[j]!=0||pf[j]==0);
   auto update_intercept=[&](){if(!intercept)return 0.0;double nom=0,den=0;for(int i=0;i<n;++i){nom+=mu[i]*r[i];den+=mu[i];}
    if(den<=0)return 0.0;double delta=nom/den;if(delta!=0){auto ts=Clock2::now();for(int i=0;i<n;++i)r[i]-=delta;tm_res+=milliseconds(ts);b0+=delta;}return std::abs(delta);};
   auto update=[&](int j){if(curvature[j]<=std::numeric_limits<double>::epsilon())return 0.0;
   auto ts=Clock2::now();double g=0;const auto& rows=support(j);++visits;coordinate_rows+=rows.size();for(int i:rows)g+=mu[i]*r[i];g*=inv;tm_grad+=milliseconds(ts);
    // Once a coordinate enters the exact active-set update loop, retain its
    // support deterministically across sweeps, IRLS iterations and lambdas.
    promote(j,rows);
    double candidate=st(curvature[j]*beta[j]+g,lam*pf[j])/curvature[j];double delta=candidate-beta[j];
    if(delta!=0){ts=Clock2::now();for(int i:rows)r[i]-=delta;tm_res+=milliseconds(ts);beta[j]=candidate;if(candidate!=0)active[j]=1;}
    return std::abs(delta);};
   // Do not perform an O(n*p) univariate support scan on the first
   // sweep. Begin with current active / unpenalized coefficients only;
   // one simultaneous cumulative-gradient KKT scan screens inactive terms.
   double change=update_intercept();
   for(int j=0;j<p;++j)if(!active_set||active[j])change=std::max(change,update(j));
   ++total_sweeps;
   bool inner_ok=false;
   for(int sweep=1;sweep<maxit_cd;++sweep){change=update_intercept();for(int j=0;j<p;++j)if(!active_set||active[j])change=std::max(change,update(j));++total_sweeps;
    if(change<=tol*(1+std::abs(b0))){auto ts=Clock2::now();bool added=false;
     if(active_set){auto tk=Clock2::now();++inactive_scans;++screening_scans;for(int i=0;i<n;++i)wr[i]=mu[i]*r[i];NumericVector gradient=full_sums(ones,wr)["gradient"];tm_inactive+=milliseconds(tk);
      for(int j=0;j<p;++j)if(eligible[j]&&!active[j]&&std::abs(gradient[j])>lam*pf[j]+tol){active[j]=1;added=true;}}
     tm_checks+=milliseconds(ts);if(!added){inner_ok=true;break;}}
    if(sweep%100==0)checkUserInterrupt();
   }
   all_inner=all_inner&&inner_ok;
   auto tl=Clock2::now();double proposed=likelihood(lam,false);tm_likelihood+=milliseconds(tl);
   if(!std::isfinite(proposed)||proposed>current+1e-12*(1+std::abs(current))){auto candidate=beta;double candidate0=b0;bool accepted=false;
    for(int ls=1;ls<=30;++ls){double fraction=std::ldexp(1.0,-ls);for(int j=0;j<p;++j)beta[j]=old[j]+fraction*(candidate[j]-old[j]);b0=old0+fraction*(candidate0-old0);
     tl=Clock2::now();proposed=likelihood(lam,false);tm_likelihood+=milliseconds(tl);if(std::isfinite(proposed)&&proposed<=current){accepted=true;break;}}
    if(!accepted){beta=old;b0=old0;proposed=current;}}
   tl=Clock2::now();double final_obj=likelihood(lam,true);tm_likelihood+=milliseconds(tl);
   // Full KKT scan covers EVERY coefficient, including screened-out ones.
   // A strong-rule mistake must be corrected by adding the violating
   // coordinate and repeating IRLS; never report convergence otherwise.
   tic=Clock2::now();
   std::vector<double> full_signal(n);for(int i=0;i<n;++i)full_signal[i]=mu[i]-y[i];
   NumericVector full_gradient=full_sums(ones,full_signal)["gradient"];
   ++strong_rule_scans[kk];
   final_kkt=0;
   for(int j=0;j<p;++j){
     double g=full_gradient[j], bound=lam*pf[j],violation;
     if(beta[j]>1e-12)violation=std::abs(g+bound);
     else if(beta[j]<-1e-12)violation=std::abs(g-bound);
     else violation=std::max(0.0,std::abs(g)-bound);
     final_kkt=std::max(final_kkt,violation);
     if(!eligible[j] && violation>std::max(10*tol,1e-8)){
       eligible[j]=1;++reentries;
     }
   }
   if(intercept){double g0=0;for(double z:full_signal)g0+=z;
     final_kkt=std::max(final_kkt,std::abs(g0*inv));}
   tm_final_kkt+=milliseconds(tic);tm_checks+=milliseconds(tic);
   if(reentries>0)strong_candidates[kk]=count_eligible();

   double delta_param=std::abs(b0-old0);for(int j=0;j<p;++j)delta_param=std::max(delta_param,std::abs(beta[j]-old[j]));
   double delta_obj=std::abs(final_obj-current)/(1+std::abs(current));iterations[kk]=outer+1;
   if(inner_ok&&delta_obj<=tol&&delta_param<=std::sqrt(tol)*(1+std::abs(b0))&&final_kkt<=std::max(10*tol,1e-8)){outer_ok=true;break;}
  }
  auto tl=Clock2::now();double obj=likelihood(lam,true),pen=0,loss=0;tm_likelihood+=milliseconds(tl);int active=0;
  for(int j=0;j<p;++j){B(j,kk)=beta[j];pen+=lam*pf[j]*std::abs(beta[j]);if(beta[j]!=0)++active;}
  for(int i=0;i<n;++i)loss+=mu[i]-y[i]*eta[i];
  a0[kk]=b0;objective_path[kk]=obj;loss_path[kk]=loss*inv;pen_path[kk]=pen;kkt_path[kk]=final_kkt;
  sweeps[kk]=total_sweeps;active_size[kk]=active;converged[kk]=outer_ok;cd_converged[kk]=all_inner;
  timing(kk,0)=tm_irls;timing(kk,1)=tm_curv;timing(kk,2)=tm_grad;timing(kk,3)=tm_res;timing(kk,4)=tm_checks;timing(kk,5)=milliseconds(full_start);timing(kk,6)=support_ms-support_ms_begin;timing(kk,7)=tm_inactive;timing(kk,8)=tm_final_kkt;timing(kk,9)=tm_likelihood;
  counts(kk,0)=double(support_builds-builds_begin);counts(kk,1)=double(support_hits-hits_begin);counts(kk,2)=double(support_rows-rows_begin);counts(kk,3)=double(visits);counts(kk,4)=double(coordinate_rows);counts(kk,5)=double(inactive_scans);
  counts(kk,6)=static_cast<double>(peak_cached_entries);counts(kk,7)=0;counts(kk,8)=static_cast<double>(screening_scans);
  strong_kkt_reentries[kk]=reentries;
  // Retain gradient only for a KKT-converged preceding lambda. The next
  // strong screen is disabled if this lambda has failed to converge.
  if(outer_ok){
   std::vector<double> gsignal(n);for(int i=0;i<n;++i)gsignal[i]=mu[i]-y[i];
   NumericVector gg=full_sums(ones,gsignal)["gradient"];
   for(int j=0;j<p;++j)preceding_gradient[j]=gg[j];
  }
 }
 return List::create(_["beta"]=B,_["a0"]=a0,_["lambda"]=lambda,_["loss"]=loss_path,_["penalty"]=pen_path,_["objective"]=objective_path,_["max_kkt"]=kkt_path,
  _["irls_iterations"]=iterations,_["cd_sweeps"]=sweeps,_["active_size"]=active_size,_["converged"]=converged,_["cd_converged"]=cd_converged,_["timing_ms"]=timing,_["profile_counts"]=counts,_["penalty_factor_scale"]=scale,
  _["strong_candidates"]=strong_candidates,_["strong_excluded"]=strong_excluded,
  _["strong_kkt_reentries"]=strong_kkt_reentries,_["strong_rule_scans"]=strong_rule_scans);
}
