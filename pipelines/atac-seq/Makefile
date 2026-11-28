# =============================================================================
#  Makefile —— 常用操作入口
#
#    make            显示帮助
#    make status     参考数据都指到哪了
#    make doctor     软件体检
#    make deploy     只补装缺失的软件
#    make test       冒烟测试（不需要生信工具）
#    make lint       语法检查
# =============================================================================
SHELL  := /bin/bash
PIPE   := bin/atac-seq
DEPLOY := deploy/deploy.sh

.PHONY: help doctor status deploy test lint clean-results

help:
	@echo "atac-seq"
	@echo ""
	@echo "  make doctor        软件体检（装了什么、还缺什么）"
	@echo "  make status        参考数据路径与大小"
	@echo "  make deploy        只补装缺失的必需工具"
	@echo "  make test          冒烟测试（不需要生信工具）"
	@echo "  make lint          语法检查"
	@echo "  make clean-results 删除 results/（不动 reference/）"
	@echo ""
	@echo "跑数据："
	@echo "  $(PIPE) -1 R1.fastq.gz -2 R2.fastq.gz -t 20"

doctor:
	@bash $(DEPLOY) doctor

status:
	@bash $(DEPLOY) status

deploy:
	@bash $(DEPLOY) env

test:
	@bash tests/smoke_test.sh

lint:
	@for f in bin/atac-seq src/*.sh deploy/*.sh tests/*.sh; do \
		printf '%-34s' "$$f"; bash -n "$$f" && echo OK || exit 1; \
	done

clean-results:
	@rm -rf results
	@echo "已删除 results/"
