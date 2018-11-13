#/bin/sh

download=0

download_extract()
{
  if which curl
  then
    echo Downloading scp.tar.xz with curl
    if curl -L -O "https://vle-project.org/pub/baryonyx/$1"
    then
      download=1
    fi
  else
    if which wget
    then
      echo Downloading scp.tar.xz with wget 
      if wget "https://vle-project.org/pub/baryonyx/$1"
      then
        download=1
      fi
    fi
  fi

  if test "$download" -eq "1"
  then
    echo Extracting scp.tar.xz
    if tar Jxf "$1" -C .
    then
      echo Success
    fi
    rm "$1"
  fi
}

download_extract scp.tar.xz
download_extract csplib022.tar.xz
download_extract spp.tar.xz
download_extract spp-wcsp2.tar.xz
download_extract telebus.tar.xz
