# =============================================================================
#  Makefile —— 常用操作入口
#
#    make            显示帮助
#    make status     查看依赖安装情况
#    make deploy     装齐依赖（等价于 bash deploy/deploy.sh all）
#    make test       冒烟测试
#    make lint       所有脚本语法检查
# =============================================================================
SHELL := /bin/bash
PIPE  := bin/germline-wes-wgs
DEPLOY := deploy/deploy.sh

.PHONY: help doctor status deploy test lint clean-results env ref bundle vep

help:
	@echo "germline-wes-wgs"
	@echo ""
	@echo "  make doctor        软件体检（已装了什么、还缺什么）"
	@echo "  make status        查看参考数据安装情况"
	@echo "  make env           只补装缺失的软件"
	@echo "  make ref           下载参考基因组 + 索引"
	@echo "  make bundle        下载 GATK 已知位点"
	@echo "  make vep           下载 VEP cache"
	@echo "  make deploy        全部安装（ANNOVAR 需手动）"
	@echo "  make test          冒烟测试（不需要生信工具）"
	@echo "  make lint          语法检查"
	@echo "  make clean-results 删除 results/（不动 reference/）"
	@echo ""
	@echo "跑数据："
	@echo "  $(PIPE) -1 R1.fq.gz -2 R2.fq.gz -c bcftools -a vep -t 16"

doctor:
	@bash $(DEPLOY) doctor

status:
	@bash $(DEPLOY) status

env:
	@bash $(DEPLOY) env

ref:
	@bash $(DEPLOY) ref

bundle:
	@bash $(DEPLOY) bundle

vep:
	@bash $(DEPLOY) vep

deploy:
	@bash $(DEPLOY) all

test:
	@bash tests/smoke_test.sh

lint:
	@for f in bin/germline-wes-wgs src/*.sh deploy/*.sh tests/*.sh; do \
		printf '%-34s' "$$f"; bash -n "$$f" && echo OK || exit 1; \
	done

clean-results:
	@rm -rf results
	@echo "已删除 results/"
