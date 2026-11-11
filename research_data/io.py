"""严格读取资料，保留来源摘要及有限数值口径。"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

def finite(x):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x):raise ValueError('需要有限数值')
    return x

def positive(x):
    if finite(x)<=0:raise ValueError('需要正数')
    return x

def integer(x,minimum=0):
    if type(x) is not int or x<minimum:raise ValueError('整数超出范围')
    return x

def identity(x):
    if not isinstance(x,str) or not x.strip():raise ValueError('身份不能为空')
    return x

def unique(rows,key):
    seen=set()
    for r in rows:
        x=identity(r[key])
        if x in seen:raise ValueError('身份重复')
        seen.add(x)
    return seen

def rational(x):
    if isinstance(x,bool) or not isinstance(x,(str,int)):raise ValueError('有理数使用整数或分数字符串')
    try:return Fraction(x)
    except (ValueError,ZeroDivisionError):raise ValueError('无效有理数') from None

def load(path):
    raw=Path(path).read_bytes()
    value=json.loads(raw.decode('utf8'),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('不接受非有限JSON常量')))
    if not isinstance(value,dict):raise ValueError('顶层需要对象')
    return value,hashlib.sha256(raw).hexdigest()
