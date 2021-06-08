#pragma once
#include <iostream>
#include <string>
using namespace std;

class Error_throw_pause 
{
    public:
    string msg;
    Error_throw_pause(string s):msg(s){};
};

class Error_throw_ignore
{
    public:
    string msg;
    Error_throw_ignore(string s):msg(s){};
};