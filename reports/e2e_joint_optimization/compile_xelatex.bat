@echo off
setlocal
set "XELATEX=D:\latex\101\texlive\2023\bin\windows\xelatex.exe"
"%XELATEX%" -interaction=nonstopmode -halt-on-error main.tex
"%XELATEX%" -interaction=nonstopmode -halt-on-error main.tex
