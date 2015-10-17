#ifndef PrematOrder_h
#define PrematOrder_h

//! Define 2 different orders on the indices of the matrices

#include <utility>
using namespace std;
template<MType T> struct cppos
{
  typedef pair<int,int> pos;
  inline bool operator()(const pos pos1,const pos pos2) const
  {
    return pos1<pos2;
  }
};
template<> struct cppos<LineWise>
{
  //type NCformat for SuperLU (SLU_NC).
  typedef pair<int,int> pos;
  inline bool operator()(const pos pos1,const pos pos2) const
  {
    return pos1.second<pos2.second||
      (pos1.second==pos2.second && 
       pos1.first<pos2.first);
  }
};
template<MType T> struct indices
{
  typedef pair<int,int> pos;
  inline int first(const pos& p)const {return p.first;}
  inline int second(const pos& p)const {return p.second;}
};
template<> struct indices<LineWise>
{
  typedef pair<int,int> pos;
  inline int first(const pos& p)const {return p.second;}
  inline int second(const pos& p)const {return p.first;}
};
#endif
