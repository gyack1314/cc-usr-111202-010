"""检查中心邻接关系并汇总计数和期望占比。"""
from .io import unique,rational,integer

def summarize(v):
    centers=v['centers'];known=unique(centers,'id');total=0;expect=0
    if not known:raise ValueError('缺少中心')
    for r in centers:
        e=rational(r['expected'])
        if e<=0:raise ValueError('期望必须正')
        expect+=e;total+=integer(r['observed'])
    edges=set()
    for edge in v['edges']:
        if len(edge)!=2 or any(x not in known for x in edge) or edge[0]==edge[1]:raise ValueError('边引用错误')
        edges.add(tuple(sorted(edge)))
    return {'centers':len(known),'undirected_edges':len(edges),'observed_total':total,'expected_total':str(expect),'expected_shares':{r['id']:str(rational(r['expected'])/expect) for r in sorted(centers,key=lambda x:x['id'])}}
