# Python 的 repr() 与双下划线方法

## 1. repr() 是什么函数

`repr()` 是 Python 的内置函数，用于取得一个对象的字符串表示，主要方便开发者查看和调试。它的返回值是字符串。

```python
repr(123)       # 返回字符串 '123'
repr([1, 2])    # 返回字符串 '[1, 2]'
```

`str()` 和 `repr()` 都返回字符串，但侧重点不同：`str()` 偏向给用户阅读，`repr()` 偏向展示对象的明确表示。

```python
s = "hello\nworld"

print(str(s))
# hello
# world

print(repr(s))
# 'hello\nworld'
```

这里 `repr(s)` 的结果保留了表示字符串的引号，并把换行显示成 `\n`，便于看清字符串实际包含什么。

对于许多内置类型，`repr()` 的结果看起来像可以重新创建对象的 Python 表达式，但并非所有对象都保证这一点。例如，DataFrame 的表示通常是供人查看的表格文本。

## 2. Python 中双下划线的作用

### 前后都有双下划线：特殊方法或属性

名称前后都有双下划线，例如 `__repr__`，通常表示 Python 预先规定的特殊方法或特殊属性，也常叫 *dunder*（double underscore）。特殊方法让你定义一个类如何响应 Python 的标准操作。

| 使用的操作 | 对应的特殊方法 | 作用 |
|---|---|---|
| `repr(obj)` | `__repr__` | 返回对象的调试字符串 |
| `str(obj)` | `__str__` | 返回对象的可读字符串 |
| `len(obj)` | `__len__` | 返回对象的长度 |
| `obj[key]` | `__getitem__` | 获取元素 |
| `a + b` | `__add__` | 定义加法行为，完整运算规则还可能涉及其他方法 |

例如，定义 `__len__` 后，自己的类就可以支持 `len()`：

```python
class Bag:
    def __len__(self):
        return 3


bag = Bag()
len(bag)  # 3
```

不是给任意方法加双下划线就会产生特殊功能，需要使用 Python 已经规定的名称。双下划线名称也不一定是方法，例如 `__name__` 是特殊属性或模块变量。

### 只有开头有双下划线：名称改写

类中只有开头有双下划线、末尾没有双下划线的名称，例如 `__value`，会触发名称改写，主要用来避免子类中的同名属性冲突。

```python
class Example:
    def __init__(self):
        self.__value = 10


obj = Example()
print(obj._Example__value)  # 10
```

这不是严格的访问权限限制。`__value` 与前后都有双下划线的 `__repr__` 属于不同机制。

## 3. repr() 和 __repr__() 的关系

`repr()` 是内置函数，`__repr__()` 是类中定义的特殊方法，两者不是同一个函数。当调用 `repr(obj)` 时，Python 会使用对象所属类的 `__repr__` 方法生成字符串；该方法必须返回字符串。

```python
class Person:
    def __repr__(self):
        return "Person(name='Shawn')"


p = Person()
print(repr(p))  # Person(name='Shawn')
```

在普通 Python 交互环境中直接输入 `p`，也会通过 `repr()` 展示对象。某些环境还支持其他显示方式，例如 Jupyter 中的 HTML 表格。

## 4. 仓库中曾使用的写法

`common/frame.py` 中的 `TimeSeriesDataFrame` 继承自 `pd.DataFrame`。目前已删除自定义的 `__repr__` 和 `_repr_html_`，直接继承 pandas 的显示方法。下面保留原来的写法，用来说明调用关系：

```python
def __repr__(self):
    """Use pandas' base formatter without subclass reconstruction."""
    return repr(pd.DataFrame(self, copy=False))
```

可以把它拆开理解：

```python
def __repr__(self):
    plain_df = pd.DataFrame(self, copy=False)
    return repr(plain_df)
```

先把当前对象包装成普通 DataFrame，再调用内置函数 `repr()`，使用普通 DataFrame 的显示方法。`copy=False` 表示尽量避免复制底层数据。

```text
repr(TimeSeriesDataFrame 对象)
    → TimeSeriesDataFrame.__repr__
        → 转成普通 pd.DataFrame
        → repr(普通 DataFrame)
            → pd.DataFrame.__repr__
                → 返回表格字符串
```

这里不会再次调用自己的 `__repr__`，因为传给内层 `repr()` 的对象已经是普通 DataFrame。

原来这层包装的意图是避免显示过程中的临时表格操作触发子类构造，并没有定制显示内容。是否需要这种处理，应根据实际问题判断，不能仅凭可能触发子类构造就认定必须重写。

如果改成下面这样，才会不断调用自己，最终触发递归错误：

```python
def __repr__(self):
    return repr(self)  # 又调用当前对象的 __repr__
```

原来相邻的 `_repr_html_` 则用于 Jupyter 等环境中的 HTML 显示。注意它前后各只有一个下划线，属于这些环境采用的显示约定，不是 Python 内置的双下划线方法；删除子类中的实现后，仍可继承 pandas 的同名方法。
