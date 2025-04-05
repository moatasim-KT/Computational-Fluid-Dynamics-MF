Traceback (most recent call last):
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/main.py", line 61, in main
    stmt_list = parse.parse(buf if buf[-1] == '\n' else buf + '\n')
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/parse.py", line 849, in parse
    p = parser.parse(
        ^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/ply/yacc.py", line 331, in parse
    return self.parseopt(input, lexer, debug, tracking, tokenfunc)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/ply/yacc.py", line 909, in parseopt
    tok = call_errorfunc(self.errorfunc, errtoken, self)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/ply/yacc.py", line 192, in call_errorfunc
    r = errorfunc(token)
        ^^^^^^^^^^^^^^^^
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/parse.py", line 836, in p_error
    raise_exception(SyntaxError,
  File "/usr/local/python/3.12.1/lib/python3.12/site-packages/smop/lexer.py", line 335, in raise_exception
    raise error_type(message, (options.filename,
  File "./FiniteDifferenceMethods/exercise/laplace_ex2.m", line 30
    % Kronecker Delta product
    ^
SyntaxError: Unexpected ";

" (parser)
Errors: 1
