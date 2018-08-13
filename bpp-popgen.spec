%define _prefix /usr

URL: https://github.com/BioPP/bpp-seq

Name: bpp-popgen
Version: 2.4.1
Release: 1%{?dist}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: %{name}-%{version}.tar.gz
Summary: Bio++ Population Genetics library
Group: Development/Libraries/C and C++
Requires: bpp-core = %{version}
Requires: bpp-seq = %{version}

BuildRoot: %{_builddir}/%{name}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: libbpp-core4 = %{version}
BuildRequires: libbpp-core-devel = %{version}
BuildRequires: libbpp-seq12 = %{version}
BuildRequires: libbpp-seq-devel = %{version}

AutoReq: yes
AutoProv: yes

%description
This library contains utilitary and classes for population genetics analysis. 
It is part of the Bio++ project.

%package -n libbpp-popgen8
Summary: Bio++ Population Genetics library
Group: Development/Libraries/C and C++

%description -n libbpp-popgen8
This library contains utilitary and classes for population genetics and molecular evolution analysis.
It is part of the Bio++ project.

%package -n libbpp-popgen-devel
Summary: Libraries, includes to develop applications with %{_basename}
Group: Development/Libraries/C and C++
Requires: libbpp-popgen8 = %{version}
Requires: libbpp-seq12 = %{version}
Requires: libbpp-seq-devel = %{version}
Requires: libbpp-core4 = %{version}
Requires: libbpp-core-devel = %{version}

%description -n libbpp-popgen-devel
The libbpp-popgen-devel package contains the header files and static libraries for
building applications which use %{name}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DBUILD_TESTING=OFF"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -n libbpp-popgen8 -p /sbin/ldconfig

%postun -n libbpp-popgen8 -p /sbin/ldconfig

%files -n libbpp-popgen8
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so.*

%files -n libbpp-popgen-devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%dir %{_prefix}/%{_lib}/cmake/
%dir %{_prefix}/%{_lib}/cmake/bpp-popgen
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/%{_lib}/cmake/bpp-popgen/bpp-popgen*.cmake
%{_prefix}/include/*

%changelog
* Mon Aug 13 2018 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.4.1-1
- Compatibility update gcc8.
* Mon Mar 12 2018 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.4.0-1
- Remove deprecated dynamic exception specifications
- Increased interface number
* Tue Jun 06 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.1-1
- Increased interface number
* Wed May 10 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.0-1
- Several bugs fixed and performance improvements
- Improved PolymorphismSequenceContainer
- Upgrade to C++11
* Fri Sep 19 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.2.0-1
- Statistical funcion names rationalized.
* Thu Mar 07 2013 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.1.0-1
- Bug fixed and warnings removed.
* Thu Feb 09 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.3-1
- Recompilation for dependencies.
* Thu Jun 09 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.2-1
- New Fst calculations + bugs fixed.
* Mon Feb 28 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.1-1
* Mon Feb 07 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.0-1
* Thu Mar 25 2010 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.5.0-1
* Wed Jun 10 2009 Julien Dutheil <jdutheil@birc.au.dk> 1.4.0-1
* Thu Dec 11 2008 Julien Dutheil <jdutheil@birc.au.dk> 1.3.1-1
* Mon Jul 21 2008 Julien Dutheil <jdutheil@birc.au.dk> 1.3.0-1
* Fri Jan 18 2008 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.2.0-1
* Fri Jul 06 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.1.1-1
- For compatibility. No more dependency for Bpp-Phyl.
* Fri Jan 19 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.1.0-2
- Build 2 for compatibility.
* Mon Aug 28 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.1.0-1
- Now requires Bpp-Phyl too!
* Tue Apr 18 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.0.0-2
- Build 2 for compatibility with other libs. Added STL dependency.
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.0.0-1
- First draft of the spec file.

