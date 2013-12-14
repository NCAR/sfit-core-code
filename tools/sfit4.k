#!/bin/ksh

print ' This script :' $0
#print $#

ver='release'
specasc=0
convert=0
clean=0
quit=0
newbin=
h2bin=0
verbos=0
nodtl=0

base=/Users/jamesw/FDP/sfit/400
#base=/data/tools/400
bin=~/bin

#base=/project/ya4/sfitsrc
#bin=/project/ya4/bin

while getopts ":b:sv:chlqzVY" opt; do
    case $opt in
        b  ) newbin=$OPTARG ;;
        s  ) specasc=1 ;;
        v  ) ver=$OPTARG ;;
        l  ) clean=1 ;;
        z  ) clean=2 ;;
        q  ) quit=1 ;;
        c  ) convert=1 ;;
        V  ) verbos=1 ;;
        Y  ) nodtl=1 ;;
        h  ) h2bin=1 ;;
        \? ) print 'usage: sfit4.k -s 0 -v 3.81 -b ~/bin -c -h -q'
             print '    -b : override binary directory'
             print '    -s : make t15asc [0|1]'
             print '    -v : version [v003.81 | v4 | test]'
             print '    -c : convert old binput'
             print '    -V : verbose output to screen'
             print '    -Y : no output to detail file'
             print '    -l : clean old sfit4 output and quit'
             print '    -h : run hbin'
             print '    -q : quit after pre-processor (hbin, prepspec, convert)'
             print '    default: do not make t15asc, do not convert binput, use version 3.81'
             return 1 ;;
    esac
done

case $ver in
   release )
       rdrv=../../src/sfit4
       pspc=../../src/pspec
       hbin=../../src/hbin
       ;;
   git )
       rdrv=${base}/sfit-core/src/sfit4
       pspc=${base}/sfit-core/src/pspec
       hbin=${base}/sfit-core/src/hbin
       ;;
   work )
       rdrv=${base}/src/delta/work/sfit4
       pspc=${base}/src/delta/work/pspec
       hbin=${base}/src/delta/work/hbin
       ;;
   mp )
       rdrv=${base}/src/mp/sfit/sfit4
       pspc=${base}/src/mp/sfit/pspec
       hbin=${base}/src/mp/sfit/hbin
       ;;
   0.9.4.1 )
       rdrv=${base}/src/delta/0.9.4.1/sfit4
       pspc=${base}/src/delta/0.9.4.1/pspec
       hbin=${base}/src/delta/0.9.4.1/hbin
       ;;
   d0.9.2 )
       rdrv=${base}/src/delta/d0.9.2/sfit4
       cvrt=${base}/src/delta/d0.9.2/convert_rdrv_394
       pspc=${base}/src/delta/d0.9.2/pspec
       hbin=${base}/src/delta/d0.9.2/hbin
       ;;
   d0.9.1 )
       rdrv=${base}/src/delta/d0.9.1/sfit4
       cvrt=${base}/src/delta/d0.9.1/convert_rdrv_394
       pspc=${base}/src/delta/d0.9.1/pspec
       hbin=${base}/src/delta/d0.9.1/hbin
       ;;
   comp )
       rdrv=${base}/src/gamma/comp/sfit4
       cvrt=${base}/src/gamma/comp/convert_rdrv_394
       pspc=${base}/src/gamma/comp/pspec
       hbin=${base}/src/gamma/comp/hbin
       ;;
   p1.7 )
       rdrv=${base}/src/gamma/p1.7/sfit4
       cvrt=${base}/src/gamma/p1.7/convert_rdrv_394
       pspc=${base}/src/gamma/p1.7/pspec
       hbin=${base}/src/gamma/p1.7/hbin
       ;;
   p1.6 )
       rdrv=${base}/src/gamma/p1.6/sfit4
       cvrt=${base}/src/gamma/p1.6/convert_rdrv_394
       pspc=${base}/src/gamma/p1.6/pspec
       hbin=${base}/src/gamma/p1.6/hbin
       ;;
   p1.5 )
       rdrv=${base}/src/gamma/p1.5/sfit4
       cvrt=${base}/src/gamma/p1.5/convert_rdrv_394
       pspc=${base}/src/gamma/p1.5/pspec
       hbin=${base}/src/gamma/p1.5/hbin
       ;;
   g1.2 )
       rdrv=${base}/src/gamma/g1.2/sfit4
       cvrt=${base}/src/gamma/g1.2/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/gamma/g1.2/hbin
      ;;
   m5.0 )
       rdrv=${base}/src/ews/m5.0/sfit4
       cvrt=${base}/src/ews/m5.0/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/ews/m5.0/hbin
      ;;
   m4.0 )
       rdrv=${base}/src/ews/m4.0/sfit4
       cvrt=${base}/src/ews/m4.0/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/ews/m4.0/hbin
      ;;
   m3.91 )
       rdrv=${base}/src/ews/m3.91/sfit4
       cvrt=${base}/src/ews/m3.91/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/ews/m3.91/hbin
      ;;
   3.992 )
       rdrv=${base}/src/b3.992/sfit4
       cvrt=${base}/src/b3.992/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/b3.992/hbin
      ;;
   3.999p )
       rdrv=${base}/src/b3.999p/sfit4
       cvrt=${base}/src/b3.999p/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/b3.999p/hbin
       ;;
   3.999 )
       rdrv=${base}/src/beta/b3.999/sfit4
       cvrt=${base}/src/beta/b3.999/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/beta/b3.999/hbin
      ;;
   3.81 )
       rdrv=${base}/src/v003.81/sfit4
       cvrt=${base}/src/v003.81/convert_rdrv_394
       pspc=${base}/v003.81/prepspec
       ;;
   3.91 )
       rdrv=${base}/src/alpha/v003.91/sfit4
       cvrt=${base}/src/alpha/v003.91/convert_rdrv_394
       pspc=${base}/preps/prepspec
      ;;
   3.92 )
       rdrv=${base}/src/v003.92/sfit4
       cvrt=${base}/src/v003.92/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/v003.92/hbin
      ;;
   3.93 )
       rdrv=${base}/src/v003.93/sfit4
       cvrt=${base}/src/v003.93/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/v003.93/hbin
      ;;

   3.94 )
       rdrv=${base}/src/v003.94/sfit4
       cvrt=${base}/src/v003.94/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/v003.94/hbin
      ;;

   3.991 )
       rdrv=${base}/src/b3.991/sfit4
       cvrt=${base}/src/b3.991/convert_rdrv_394
       pspc=${base}/preps/prepspec
       hbin=${base}/src/b3.991/hbin
      ;;

   v4 )
       rdrv=${base}/src/v004.00/sfit4
       cvrt=${base}/src/v004.00/convert_rdrv_394
       pspc=${base}/v004.00/prepspec
      ;;
   test )
      rdrv=""
      cvrt=""
      pspc=""
      ;;
esac


# if bin is set then override 'version' and use this bin directory
if [ ! -z ${newbin} ]
then
   rdrv=${newbin}/sfit4
   cvrt=${newbin}/convert_rdrv_394
   pspc=${newbin}/prepspec
   hbin=${newbin}/hbin
fi

# convert old binput then stop
if ((${convert}==1))
then
   print "    Running convert_rdrv_394 version ${ver}...${cvrt}"
   (${cvrt})
   cvOK=$?
   if [ $cvOK -ne 0 ]
   then
      print '${cvrt} error returned : ' $spOK ' exiting 2'
      exit 2
   fi
fi

# run prepspec
if ((${specasc}==1))
then
   print "    Running pspec version ${ver}...${pspc} ... saving output to pspec.dtl"
   if ((${nodtl}==1))
   then
      (${pspc})
   else
      (${pspc}) > pspec.dtl
   fi
   spOK=$?
   if [ $spOK -ne 0 ]
   then
      print 'prepspec error returned : ' $spOK ' exiting 3'
      exit 3
   fi
   if ((${verbos}==1 && ${nodtl}==0))
   then
      cat pspec.dtl
   fi
fi

# run hbin
if ((${h2bin}==1))
then
   print "    Running hbin version ${ver}...${hbin} ... saving output to hbin.dtl"
   if ((${nodtl}==1))
   then
      (${hbin})
   else
      (${hbin}) > hbin.dtl
   fi
   spOK=$?
   if [ $spOK -ne 0 ]
   then
      print 'hbin error returned : ' $spOK ' exiting 3'
      exit 3
   fi
   if ((${verbos}==1 && ${nodtl}==0))
   then
      cat hbin.dtl
   fi
fi

if ((${clean}>=1))
then
   rm -f A-S.out AK.out detail detail.opt  K.out newratio. Parm.out partkrnl.dat pbpfile
   rm -f PRFS.out raytrace.out raytrace.pnch raytrace.sa Sa.complete SM.out solspec.dat
   rm -f SS.out statevec summary summary.st synspec.out raytrace.mix raytrace.ms raytrace.pt
   rm -f gas.* allgases.*.*  pltfits.ps K.out.Tp detail.lm detail.voigt, spc.*.*.* SA.out
   rm -f rprfs.out k.output fort.* rprfs.output kb.out *.out *.output sa.complete *.target
   rm -f aprfs.table rprfs.table shat.complete sainv.complete
fi

if ((${clean}>=2))
then
   rm -f sfit4.dtl pspec.dtl hbin.dtl
fi

# end early after one of the pre-processors
if ((${quit}==1))
then
   print 'sfit4.k...Done.'
   exit 0
fi


print "    Running SFit4 version $ver ${rdrv}..."
print ""
(${rdrv})
sfOK=$?
pwd

if [ $sfOK -ne 0 ]
then
   print 'sfit4 error returned : ' $sfOK ' exiting 4'
   exit 4
fi

#cat zsum.t summary.st >zsum

exit 0
