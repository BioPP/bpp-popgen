%define name bpp-popgen
%define version 1.3.1
%define release 1
%define _prefix /usr/local

Summary: The Bio++ PopGenLib library.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: The Bio++ Project
Source: http://kimura.univ-montp2.fr/BioPP/Repositories/sources/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: System Environment/Libraries
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
AutoReqProv: no
Requires: libstdc++6
Requires: bpp-utils >= 1.3.0
Requires: bpp-numcalc >= 1.6.0
Requires: bpp-seq >= 1.5.0

%description
This library contains utilitary and classes for population genetics analysis.
It is part of the Bio++ project.

%package devel
Summary: Libraries, includes to develop applications with %{name}.
Group: Development/Libraries
Requires: %{name} = %{version}
Requires: bpp-utils-devel >= 1.3.0
Requires: bpp-numcalc-devel = 1.6.0
Requires: bpp-seq-devel = 1.5.0

%description devel
The %{name}-devel package contains the header files and static libraries for
building applications which use %{name}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS" ./configure --prefix=%{_prefix}
make

%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS COPYING INSTALL NEWS README ChangeLog
%{_prefix}/lib/lib*.so
%{_prefix}/lib/lib*.so.*

%files devel
%defattr(-,root,root)
%doc AUTHORS COPYING INSTALL NEWS README ChangeLog
%{_prefix}/lib/lib*.a
%{_prefix}/lib/lib*.la
%{_prefix}/include/*

%changelog
* Thu Dec 11 2008 Julien Dutheil <jdutheil@daimi.au.dk>
- Version 1.3.1.
* Mon Jul 21 2008 Julien Dutheil <jdutheil@daimi.au.dk>
- Version 1.3.0.
* Fri Jan 18 2008 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.2.0.
* Fri Jul 06 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.1.1 for compatibility. No more dependency for Bpp-Phyl.
* Fri Jan 19 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.1.0 build 2 for compatibility.
* Mon Aug 28 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.1.0, now requires Bpp-Phyl too!
* Tue Apr 18 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Build 2 for compatibility with other libs. added STL dependency.
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- First draft of the spec file

