set foldcolumn=2
" set foldexpr=FoldPythonFuncdefs()
" set foldtext=FoldTextGeneric()
" setlocal foldmethod=expr
" setlocal foldexpr=MarkdownFolds()

" {{{1 
function! FoldPythonFuncdefs()
  let thisline = getline(v:lnum)
  if match(thisline, '^def .*:') >= 0
    return ">1"
  else
    return "="
  endif
endfunction

" {{{1
function! FoldSnakemake()
  let thisline = getline(v:lnum)
  if match(thisline, '^rule .*:') >= 0
    return ">1"
  elseif match(thisline, '^#--- ') >= 0
    return ">1"
  elseif match(thisline, '# =====') >= 0
    return ">0"
  else
    return "="
  endif
endfunction

" {{{1
function! FoldRmarkdown()
  let thisline = getline(v:lnum)
  if match(thisline, 'title:') >= 0
    return ">1"
  elseif match(thisline, '^```{r .*}') >= 0
    return ">1"
  else
    return "="
  endif
endfunction


" function! FoldTextGeneric()
"   let foldsize = (v:foldend-v:foldstart)
"   return getline(v:foldstart).' ('.foldsize.' Lines Folded)'
" endfunction
" set foldmethod=expr
" set foldexpr=FoldSnakemake()
