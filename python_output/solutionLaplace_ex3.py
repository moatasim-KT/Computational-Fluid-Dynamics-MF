str
Traceback (most recent call last):
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/main.py", line 66, in main
    G = resolve.resolve(stmt_list)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/resolve.py", line 54, in resolve
    u = G.node[n]["ident"]
        ^^^^^^
AttributeError: 'DiGraph' object has no attribute 'node'. Did you mean: '_node'?
Errors: 1
