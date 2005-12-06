[PIP2-like future input] Please enter:
- the context matrix,
0 4
- the bignum column (start at 0, -1 if no bignum),
-1
- the constraint matrix.
7 6
   1   2   6   0   0  -9
   1   5  -3   0   0   0
   1   2 -10   0   0  15
   1  -2   6   0   0  -3
   1  -2  -6   0   0  17
   1   0   1  -1   0   0
   1   1   0   0  -1   0

(if #[ 0 -4 3]
 (if #[ -4 0 5]
  ()
  (if #[ -1 0 3]
   (if #[ -1 0 2]
    ()
    ()
   )
   ()
  )
 )
 (if #[ 0 -2 9]
  (if #[ -6 -2 9]
   (if #[ 0 -2 3]
    (newparm 2 (div #[ 0 2 3] 6))
    (newparm 3 (div #[ 0 2 10 7] 12))
    (newparm 4 (div #[ 0 4 0 2 1] 6))
    ()
    (if #[ 0 -2 7]
     (newparm 2 (div #[ 0 4 3] 6))
     (if #[ 0 -8 6 11]
      ()
      ()
     )
     ()
    )
   )
   (if #[ -6 -2 17]
    (if #[ -3 5 0]
     (if #[ 6 -2 -3]
      ()
      (if #[ 0 -2 7]
       (newparm 2 (div #[ 0 4 3] 6))
       (if #[ 0 -8 6 11]
        ()
        ()
       )
       ()
      )
     )
     ()
    )
    ()
   )
  )
  ()
 )
)
