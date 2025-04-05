Traceback (most recent call last):
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/main.py", line 66, in main
    G = resolve.resolve(stmt_list)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/resolve.py", line 51, in resolve
    G = as_networkx(t)
        ^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/resolve.py", line 39, in as_networkx
    if u.lexpos < v.lexpos:
       ^^^^^^^^^^^^^^^^^^^
TypeError: '<' not supported between instances of 'NoneType' and 'NoneType'
Errors: 1
