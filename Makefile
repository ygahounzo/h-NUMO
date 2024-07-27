RM = /bin/rm -rf
MV = /bin/mv -f
CP = /bin/cp -f

export PLATFORMS:= hnumo hnumo-debug hnumo-brew hnumo-prof

export NUMO_DIR    := $(CURDIR)
export DEPEND_FILE := $(NUMO_DIR)/depend.mk
export CONFIG_USER := $(NUMO_DIR)/config.user
export P4EST_LOC   := $(NUMO_DIR)/p4est
export P4EST_VER   := 2.8
export P4EST_DIR   := $(NUMO_DIR)/p4est-$(P4EST_VER)

export TARGET=$(NUMO_DIR)/bin/numo3d

include $(NUMO_DIR)/config.p4est

SUBDIR := $(NUMO_DIR)/bin
SUBDIR += $(NUMO_DIR)/include
SUBDIR += $(NUMO_DIR)/libs

p4est_web=https://p4est.github.io/release/p4est-$(P4EST_VER).tar.gz

.NOTPARALLEL:

.PHONY: p4est

$(PLATFORMS): %: %-config p4est/local/lib/libp4est.a depend $(SUBDIR)
	$(MAKE) -C src $@ TARGET=$(TARGET)

$(SUBDIR):
	if [ ! -d $@ ]; then mkdir -p $@;fi

p4est/local/lib/libp4est.a:
	cd $(NUMO_DIR) && \
	curl -O $(p4est_web) && \
	tar -zxf p4est-$(P4EST_VER).tar.gz && \
	rm p4est-$(P4EST_VER).tar.gz && \
	mv $(P4EST_DIR) $(P4EST_LOC) && \
	cd $(P4EST_LOC) && \
	./configure $(P4EST_CONF) && \
	make -j && \
	make install

p4est_clean:
	cd $(NUMO_DIR)/p4est && \
	rm -rf * && \
	rm -r ../p4est

depend:
	$(NUMO_DIR)/make_depend.pl -s -l -o $(DEPEND_FILE) $(NUMO_DIR)/src

emptyrule:

clean : depend
	$(MAKE) -C src $@
	for dir in $(SUBDIR); do (echo $(RM) $$dir ; $(RM) $$dir) ; done

deepclean : clean depend p4est_clean
	$(MAKE) -C src $@
	(if [ -e $(DEPEND_FILE) ] ; then $(RM) $(DEPEND_FILE) ; fi )
