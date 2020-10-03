# echo "cleanup DStarIDAlg 00-00-01 in /afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtDStarIDAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtDStarIDAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=DStarIDAlg -version=00-00-01 -path=/afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics  $* >${cmtDStarIDAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=DStarIDAlg -version=00-00-01 -path=/afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics  $* >${cmtDStarIDAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtDStarIDAlgtempfile}
  unset cmtDStarIDAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtDStarIDAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtDStarIDAlgtempfile}
unset cmtDStarIDAlgtempfile
exit $cmtcleanupstatus

