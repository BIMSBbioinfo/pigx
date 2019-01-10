let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +1 pigx-bsseq.in
badd +1 BSseq_pipeline.py
badd +1 scripts/func_defs.py
argglobal
silent! argdel *
$argadd pigx-bsseq.in
$argadd BSseq_pipeline.py
set stal=2
edit pigx-bsseq.in
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winheight=1 winminwidth=1 winwidth=1
wincmd =
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
5,23fold
46,53fold
72,136fold
140,143fold
145,146fold
148,163fold
165,210fold
212,258fold
260,301fold
312,321fold
304,321fold
323,326fold
328,380fold
382,393fold
399,408fold
409,425fold
428,458fold
480,485fold
486,489fold
491,492fold
493,495fold
477,512fold
514,515fold
517,522fold
524,540fold
let s:l = 139 - ((2 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
139
normal! 011|
wincmd w
argglobal
2argu
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
23,40fold
42,50fold
59,140fold
142,146fold
52,155fold
171,196fold
201,238fold
244,278fold
284,329fold
335,357fold
363,375fold
381,405fold
411,423fold
426,435fold
441,448fold
451,458fold
466,489fold
491,517fold
523,539fold
545,557fold
563,576fold
578,594fold
600,616fold
618,639fold
645,660fold
52
normal! zo
52
normal! zc
let s:l = 618 - ((77 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
618
normal! 0
wincmd w
2wincmd w
wincmd =
tabedit scripts/func_defs.py
set splitbelow splitright
wincmd t
set winminheight=1 winheight=1 winminwidth=1 winwidth=1
argglobal
1argu
if bufexists('scripts/func_defs.py') | buffer scripts/func_defs.py | else | edit scripts/func_defs.py | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
3,21fold
23,24fold
26,32fold
34,38fold
40,48fold
50,58fold
60,67fold
69,78fold
80,87fold
89,96fold
98,103fold
105,110fold
112,119fold
121,128fold
130,132fold
134,147fold
149,162fold
164,165fold
168,184fold
186,187fold
189,193fold
197,202fold
205,220fold
223,226fold
229,254fold
let s:l = 255 - ((157 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
255
normal! 0
tabnext 1
set stal=1
if exists('s:wipebuf')
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToO
set winminheight=1 winminwidth=1
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
