#set allsel [atomselect 0 "all" frame 1000]
#measure sasa 1.5 $allsel

set nf [molinfo 0 get numframes]
set mols [atomselect 0 "all" frame 1000]
for { set i 1 } { $i <= $nf } { incr i } {
    $mols frame $i
    measure sasa 1.5 $allsel



