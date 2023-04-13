
One of the goals of this package is to have code which is extremely easy to read and which looks just like the equations in the paper. As such, we try to adhere to the following code-style as much as possible. In order of precendence, you should try to 1) follow this style guide, 2) follow the style of surrounding code, 3) consider https://github.com/invenia/BlueStyle (which also inspired this guide, although is not exactly the same). 

# Spaces

Spaces help improve readiability and should be used between arguments, in function definitions, between infix operators, etc... E.g.:

```julia
# Yes: 
foo(x, y) = x + y
bar(x, y)

# No: 
foo(x,y)=x+y
bar(x,y)
```

In long expressions, its ok to omit spaces for compactness but it should be only on the inner-most depths and consistent across expressions. E.g.:

```julia
# Ok, if the expression is really long:
foo(bar(x,y), baz(z,w))

# No, spaces not on inner-most depth:
foo(bar(x, y),baz(z, w))

# No, spaces at inconsistent depths:
foo(bar(x, y), baz(z,w))
```

No spaces around the `::` in function arguments or in their type parameters (since that breaks the "inner-most" depth thing):

```julia
# Yes:
foo(x::T, y::S) = ...

# No:
foo(x :: T, y :: S) = ...
```


# Multiline expressions

Try to keep expressions visible on one screen, about 100 columns max. 

For function definitions which require wrapping, put one argument per line with only one extra indentation level, and the closing parenthesis on the original indentation level:

```julia

# Yes:
function few_args(x, y)
    # ...
end

# Yes: 
function lots_of_args(
    arg = 1;
    kwarg = 2,
    kwargs...
)

# Yes: 
function only_kwargs(; # <- the ';' goes here
    kwarg = 2,
    kwargs...
)


# No:
function foo(
    arg = 1;
    kwarg = 2,
    kwargs...
    )

# No
function foo(
    x, y,
    z
```

Similarly, for funtion calls or multi-line expressions inside code, the inner code should be indented just one level with parenthesis taking up their own line, and the closing parenthesis should be unindented:

```julia
# Yes:
x = foo(
    a * b,
    (
        bar(c, d) +
        e + f
    )
)

# No:
x = foo(
        a * b, 
        (bar(c, d) +
         e + f))
```

