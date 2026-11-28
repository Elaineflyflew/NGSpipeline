# =============================================================================
#  Makefile —— 常用操作入口
#
#    make            显示帮助
#    make doctor     软件体检
#    make status     查看参考数据安装情况
#    make deploy     装齐依赖
#    make test       冒烟测试（不需要生信工具）
#    make lint       所有脚本语法检查
# =============================================================================
SHELL  := /bin/bash
PIPE   := bin/atac
DEPLOY := deploy/deploy.sh

.PHONY: help doctor status deploy test lint clean-results \
        env genome bt2 annotation blacklist homer atacqc igv

help:
	@echo "atac —— ATAC-seq 双端 FASTQ 分析流水线"
	@echo ""
	@echo "  make doctor        软件体检（已装了什么、还缺什么）"
	@echo "  make status        查看参考数据安装情况"
	@echo "  make env           只补装缺失的软件"
	@echo "  make genome        下载参考基因组 fasta"
	@echo "  make bt2           建 bowtie2 索引（约 1 小时）"
	@echo "  make annotation    下载 GENCODE GTF"
	@echo "  make blacklist     下载 ENCODE 黑名单"
	@echo "  make homer         安装 HOMER + hg38（⑫⑬ 要用）"
	@echo "  make atacqc        安装 ATACseqQC（可选）"
	@echo "  make igv           下载 IGV（可选）"
	@echo "  make deploy        全部安装"
	@echo "  make test          冒烟测试 + 单元测试"
	@echo "  make lint          语法检查"
	@echo "  make clean-results 删除 results/（不动 reference/）"
	@echo ""
	@echo "跑数据："
	@echo "  $(PIPE) -1 raw/R1.fastq.gz -2 raw/R2.fastq.gz -t 20"

doctor:
	@bash $(DEPLOY) doctor

status:
	@bash $(DEPLOY) status

env:
	@bash $(DEPLOY) env

genome:
	@bash $(DEPLOY) genome

bt2:
	@bash $(DEPLOY) bt2

annotation:
	@bash $(DEPLOY) annotation

blacklist:
	@bash $(DEPLOY) blacklist

homer:
	@bash $(DEPLOY) homer

atacqc:
	@bash $(DEPLOY) atacqc

igv:
	@bash $(DEPLOY) igv

deploy:
	@bash $(DEPLOY) all

test:
	@bash tests/smoke_test.sh
	@bash tests/run_tests.sh

lint:
	@for f in bin/atac src/*.sh deploy/*.sh config/atac.conf tests/*.sh; do \
		printf '%-34s' "$$f"; bash -n "$$f" && echo OK || exit 1; \
	done

clean-results:
	@rm -rf results
	@echo "已删除 results/"
