./fsplice Human.dat seq.fa -seq:5 -D:0 -thr_p:90

where:
" ./fsplice"        - name of the executable file
"Human.dat"         - parameter file for specified organizm.
"seq.fa"            - name of the input data file ("seq.fa")
"-seq:5"            - print spice site sequence, N nucleotides left and right.
"-D:0"              - dir = 0, search in direct chain only. (default)
					  dir = 1, search in reverse chain only.
                      dir = 2, search in both chains.
"-thr_p:90"         - set threshold fo all sites, to find more then N persent of true sites.

