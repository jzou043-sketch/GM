# Fetch the native third-party binary at the pinned official repository commit.
aa<-commandArgs(FALSE);script<-sub('^--file=','',aa[grepl('^--file=',aa)][1])
root<-dirname(dirname(normalizePath(script,winslash='/')))
if(!requireNamespace('jsonlite',quietly=TRUE))stop('Install jsonlite to read the binary manifest')
m<-jsonlite::fromJSON(file.path(root,'metadata/alphamate_binary_sources.json'))
name<-if(.Platform$OS.type=='windows')'AlphaMate.exe' else 'AlphaMate_Unix'
r<-m[m$File==name,];stopifnot(nrow(r)==1L)
target<-file.path(root,'code/original/vendor',name)
if(file.exists(target)) {
  stopifnot(unname(tools::md5sum(target))==r$MD5)
} else {
  partial<-paste0(target,'.download');if(file.exists(partial))stop('Inspect existing partial download: ',partial)
  download.file(r$URL,partial,mode='wb',method='libcurl')
  stopifnot(file.info(partial)$size==r$Bytes,unname(tools::md5sum(partial))==r$MD5)
  if(!file.rename(partial,target))stop('Could not save verified binary')
}
if(.Platform$OS.type!='windows')Sys.chmod(target,mode='0755')
cat('Verified native executable:',target,'\n')
