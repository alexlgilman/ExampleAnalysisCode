# echo "setup DStarIDAlg 00-00-01 in /afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtDStarIDAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtDStarIDAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=DStarIDAlg -version=00-00-01 -path=/afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics  -no_cleanup $* >${cmtDStarIDAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=DStarIDAlg -version=00-00-01 -path=/afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics  -no_cleanup $* >${cmtDStarIDAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtDStarIDAlgtempfile}
  unset cmtDStarIDAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtDStarIDAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtDStarIDAlgtempfile}
unset cmtDStarIDAlgtempfile
return $cmtsetupstatus

