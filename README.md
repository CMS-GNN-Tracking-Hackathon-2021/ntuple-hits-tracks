# ExtractHitsTracks package

```
cmsrel CMSSW_11_2_4
cd CMSSW_11_2_4/src
cmsenv
git clone git@github.com:bainbrid/ExtractHitsTracks.git test
scram b
cd test/ExtractHitsTracks/run
voms-proxy-init --voms cms
cmsRun test_cfg.py
```
