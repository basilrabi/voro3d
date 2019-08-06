#!/bin/bash

rm -rf .Rhistory \
       man \
       NAMESPACE \
       R/RcppExports.R \
       src/*o \
       src/RcppExports*

touch NAMESPACE
