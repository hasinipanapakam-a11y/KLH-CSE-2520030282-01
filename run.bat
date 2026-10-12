@echo off
rem Compile and run. Extra arguments: test | demo | scale
if not exist out mkdir out
javac -d out -sourcepath src src\texthack\app\Main.java
if errorlevel 1 exit /b 1
java -Xmx2g -cp out texthack.app.Main %*
