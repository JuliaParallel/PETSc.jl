# autodefined type arguments for class ------
mutable struct _n_SSL_CTX end
const SSL_CTX = Ptr{_n_SSL_CTX}

mutable struct _n_SSL end
const SSL = Ptr{_n_SSL}

mutable struct _n_hid_t end
const hid_t = Ptr{_n_hid_t}

mutable struct _n_MPI_Request end
const MPI_Request = Ptr{_n_MPI_Request}

const PetscDLHandle = Ptr{Cvoid}

#mutable struct _n_PetscClassId end
#const PetscClassId = Ptr{_n_PetscClassId}

#mutable struct _n_PetscInfoCommFlag end
#const PetscInfoCommFlag = Ptr{_n_PetscInfoCommFlag}

#mutable struct _n_PetscDLMode end
#const PetscDLMode = Ptr{_n_PetscDLMode}

#mutable struct _n_PetscFPTrap end
#const PetscFPTrap = Ptr{_n_PetscFPTrap}

mutable struct _n_PetscLogHandler end
const PetscLogHandler = Ptr{_n_PetscLogHandler}

#mutable struct _n_PetscLogState end
#const PetscLogState = Ptr{_n_PetscLogState}

#mutable struct _n_PetscLogStage end
#const PetscLogStage = Ptr{_n_PetscLogStage}

#mutable struct _n_PetscLogEvent end
#const PetscLogEvent = Ptr{_n_PetscLogEvent}

#mutable struct _n_PetscEventPerfInfo end
#const PetscEventPerfInfo = Ptr{_n_PetscEventPerfInfo}

mutable struct _n_PetscLayout end
const PetscLayout = Ptr{_n_PetscLayout}

mutable struct _n_PetscQuadrature end
const PetscQuadrature = Ptr{_n_PetscQuadrature}

#mutable struct _n_PetscDTSimplexQuadratureType end
#const PetscDTSimplexQuadratureType = Ptr{_n_PetscDTSimplexQuadratureType}

#mutable struct _n_DMPolytopeType end
#const DMPolytopeType = Ptr{_n_DMPolytopeType}

#mutable struct _n_LandauStaticData end
#const LandauStaticData = Ptr{_n_LandauStaticData}

mutable struct _n_P4estVertexMaps end
const P4estVertexMaps = Ptr{_n_P4estVertexMaps}

mutable struct _n_pointInterpolationP4est end
const pointInterpolationP4est = Ptr{_n_pointInterpolationP4est}

mutable struct _n_hCsize_t end
const hCsize_t = Ptr{_n_hCsize_t}

mutable struct _n_MPIU_Count end
const MPIU_Count = Ptr{_n_MPIU_Count}

mutable struct _n_PetscBLASInt end
const PetscBLASInt = Ptr{_n_PetscBLASInt}

#mutable struct _n_PetscCuBLASInt end
#const PetscCuBLASInt = Ptr{_n_PetscCuBLASInt}

#mutable struct _n_PetscHipBLASInt end
#const PetscHipBLASInt = Ptr{_n_PetscHipBLASInt}

#mutable struct _n_PetscBool3 end
#const PetscBool3 = Ptr{_n_PetscBool3}

# -------------------------------------------------------
"""
	PetscGlobusAuthorize(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, tokensize::Csize_t) 
Get an access token allowing PETSc applications to make Globus file transfer requests

Not Collective, only the first process in `MPI_Comm` does anything

Input Parameters:
- `comm`      - the MPI communicator
- `tokensize` - size of the token array

Output Parameter:
- `access_token` - can be used with `PetscGlobusUpLoad()` for 30 days

Level: intermediate

-seealso: `PetscGoogleDriveRefresh()`, `PetscGoogleDriveUpload()`, `PetscGlobusUpload()`

# External Links
$(_doc_external("Sys/PetscGlobusAuthorize"))
"""
function PetscGlobusAuthorize(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, tokensize::Csize_t) end

@for_petsc function PetscGlobusAuthorize(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, tokensize::Csize_t )

    @chk ccall(
               (:PetscGlobusAuthorize, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Csize_t),
               comm, access_token, tokensize,
              )


	return nothing
end 

"""
	PetscGlobusGetTransfers(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, buff::Vector{Cchar}, buffsize::Csize_t) 
Get a record of current transfers requested from Globus

Not Collective, only the first process in `MPI_Comm` does anything

Input Parameters:
- `comm`         - the MPI communicator
- `access_token` - Globus access token, if `NULL` will check in options database for -globus_access_token XXX otherwise
will call `PetscGlobusAuthorize()`.
- `buffsize`     - size of the buffer

Output Parameter:
- `buff` - location to put Globus information

Level: intermediate

-seealso: `PetscGoogleDriveRefresh()`, `PetscGoogleDriveUpload()`, `PetscGlobusUpload()`, `PetscGlobusAuthorize()`

# External Links
$(_doc_external("Sys/PetscGlobusGetTransfers"))
"""
function PetscGlobusGetTransfers(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, buff::Vector{Cchar}, buffsize::Csize_t) end

@for_petsc function PetscGlobusGetTransfers(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, buff::Vector{Cchar}, buffsize::Csize_t )

    @chk ccall(
               (:PetscGlobusGetTransfers, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               comm, access_token, buff, buffsize,
              )


	return nothing
end 

"""
	PetscGlobusUpload(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar}) 
Loads a file to Globus

Not Collective, only the first process in the `MPI_Comm` uploads the file

Input Parameters:
- `comm`         - MPI communicator
- `access_token` - obtained with `PetscGlobusAuthorize()`, pass `NULL` to use `-globus_access_token XXX` from the PETSc database
- `filename`     - file to upload

Options Database Key:
- `-globus_access_token XXX` - the Globus token

Level: intermediate

-seealso: `PetscGoogleDriveAuthorize()`, `PetscGoogleDriveRefresh()`, `PetscGlobusAuthorize()`

# External Links
$(_doc_external("Sys/PetscGlobusUpload"))
"""
function PetscGlobusUpload(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar}) end

@for_petsc function PetscGlobusUpload(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar} )

    @chk ccall(
               (:PetscGlobusUpload, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}),
               comm, access_token, filename,
              )


	return nothing
end 

"""
	PetscGoogleDriveRefresh(petsclib::PetscLibType,comm::MPI_Comm, refresh_token::Vector{Cchar}, access_token::Vector{Cchar}, tokensize::Csize_t) 
Get a new authorization token for accessing Google drive from PETSc from a refresh token

Not Collective, only the first process in the `MPI_Comm` does anything

Input Parameters:
- `comm`          - MPI communicator
- `refresh_token` - obtained with `PetscGoogleDriveAuthorize()`, if NULL PETSc will first look for one in the options data
if not found it will call `PetscGoogleDriveAuthorize()`
- `tokensize`     - size of the output string access_token

Output Parameter:
- `access_token` - token that can be passed to `PetscGoogleDriveUpload()`

Options Database Key:
- `-google_refresh_token XXX` - where XXX was obtained from `PetscGoogleDriveAuthorize()`

Level: intermediate

-seealso: `PetscGoogleDriveAuthorize()`, `PetscGoogleDriveUpload()`

# External Links
$(_doc_external("Sys/PetscGoogleDriveRefresh"))
"""
function PetscGoogleDriveRefresh(petsclib::PetscLibType, comm::MPI_Comm, refresh_token::Vector{Cchar}, access_token::Vector{Cchar}, tokensize::Csize_t) end

@for_petsc function PetscGoogleDriveRefresh(petsclib::$UnionPetscLib, comm::MPI_Comm, refresh_token::Vector{Cchar}, access_token::Vector{Cchar}, tokensize::Csize_t )

    @chk ccall(
               (:PetscGoogleDriveRefresh, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               comm, refresh_token, access_token, tokensize,
              )


	return nothing
end 

"""
	PetscGoogleDriveUpload(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar}) 
Loads a file to the Google Drive

Not Collective, only the first process in the `MPI_Comm` uploads the file

Input Parameters:
- `comm`         - MPI communicator
- `access_token` - obtained with PetscGoogleDriveRefresh(), pass `NULL` to have PETSc generate one
- `filename`     - file to upload; if you upload multiple times it will have different names each time on Google Drive

Options Database Key:
- `-google_refresh_token XXX` - pass the access token for the operation

-seealso: `PetscGoogleDriveAuthorize()`, `PetscGoogleDriveRefresh()`

# External Links
$(_doc_external("Sys/PetscGoogleDriveUpload"))
"""
function PetscGoogleDriveUpload(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar}) end

@for_petsc function PetscGoogleDriveUpload(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar} )

    @chk ccall(
               (:PetscGoogleDriveUpload, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}),
               comm, access_token, filename,
              )


	return nothing
end 

"""
	PetscGoogleDriveAuthorize(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, refresh_token::Vector{Cchar}, tokensize::Csize_t) 
Get authorization and refresh token for accessing Google drive from PETSc

Not Collective, only the first process in `MPI_Comm` does anything

Input Parameters:
- `comm`      - the MPI communicator
- `tokensize` - size of the token arrays

Output Parameters:
- `access_token`  - can be used with `PetscGoogleDriveUpload()` for this one session
- `refresh_token` - can be used for ever to obtain new access_tokens with `PetscGoogleDriveRefresh()`, guard this like a password
it gives access to your Google Drive

Level: intermediate

-seealso: `PetscGoogleDriveRefresh()`, `PetscGoogleDriveUpload()`

# External Links
$(_doc_external("Sys/PetscGoogleDriveAuthorize"))
"""
function PetscGoogleDriveAuthorize(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, refresh_token::Vector{Cchar}, tokensize::Csize_t) end

@for_petsc function PetscGoogleDriveAuthorize(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, refresh_token::Vector{Cchar}, tokensize::Csize_t )

    @chk ccall(
               (:PetscGoogleDriveAuthorize, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               comm, access_token, refresh_token, tokensize,
              )


	return nothing
end 

"""
	PetscSSLInitializeContext(petsclib::PetscLibType,octx::SSL_CTX) 
Set up an SSL context suitable for initiating HTTPS requests.

Output Parameter:
- `octx` - the SSL_CTX to be passed to `PetscHTTPSConnect90`

Level: advanced

If PETSc was ./configure -with-ssl-certificate requires the user have created a self-signed certificate with
-seealso: `PetscSSLDestroyContext()`, `PetscHTTPSConnect()`, `PetscHTTPSRequest()`

# External Links
$(_doc_external("Sys/PetscSSLInitializeContext"))
"""
function PetscSSLInitializeContext(petsclib::PetscLibType, octx::SSL_CTX) end

@for_petsc function PetscSSLInitializeContext(petsclib::$UnionPetscLib, octx::SSL_CTX )

    @chk ccall(
               (:PetscSSLInitializeContext, $petsc_library),
               PetscErrorCode,
               (SSL_CTX,),
               octx,
              )


	return nothing
end 

"""
	PetscSSLDestroyContext(petsclib::PetscLibType,ctx::SSL_CTX) 
frees a `SSL_CTX` obtained with `PetscSSLInitializeContext()`

Input Parameter:
- `ctx` - the `SSL_CTX`

Level: advanced

-seealso: `PetscSSLInitializeContext()`, `PetscHTTPSConnect()`

# External Links
$(_doc_external("Sys/PetscSSLDestroyContext"))
"""
function PetscSSLDestroyContext(petsclib::PetscLibType, ctx::SSL_CTX) end

@for_petsc function PetscSSLDestroyContext(petsclib::$UnionPetscLib, ctx::SSL_CTX )

    @chk ccall(
               (:PetscSSLDestroyContext, $petsc_library),
               PetscErrorCode,
               (Ptr{SSL_CTX},),
               ctx,
              )


	return nothing
end 

"""
	PetscHTTPSRequest(petsclib::PetscLibType,type::Vector{Cchar}, url::Vector{Cchar}, header::Vector{Cchar}, ctype::Vector{Cchar}, body::Vector{Cchar}, ssl::SSL, buff::Vector{Cchar}, buffsize::Csize_t) 
Send a request to an HTTPS server

Input Parameters:
- `type`     - either "POST" or "GET"
- `url`      - URL of request host/path
- `header`   - additional header information, may be `NULL`
- `ctype`    - data type of body, for example application/json
- `body`     - data to send to server
- `ssl`      - obtained with `PetscHTTPSConnect()`
- `buffsize` - size of buffer

Output Parameter:
- `buff` - everything returned from server

Level: advanced

-seealso: `PetscHTTPRequest()`, `PetscHTTPSConnect()`, `PetscSSLInitializeContext()`, `PetscSSLDestroyContext()`, `PetscPullJSONValue()`

# External Links
$(_doc_external("Sys/PetscHTTPSRequest"))
"""
function PetscHTTPSRequest(petsclib::PetscLibType, type::Vector{Cchar}, url::Vector{Cchar}, header::Vector{Cchar}, ctype::Vector{Cchar}, body::Vector{Cchar}, ssl::SSL, buff::Vector{Cchar}, buffsize::Csize_t) end

@for_petsc function PetscHTTPSRequest(petsclib::$UnionPetscLib, type::Vector{Cchar}, url::Vector{Cchar}, header::Vector{Cchar}, ctype::Vector{Cchar}, body::Vector{Cchar}, ssl::SSL, buff::Vector{Cchar}, buffsize::Csize_t )

    @chk ccall(
               (:PetscHTTPSRequest, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{SSL}, Ptr{Cchar}, Csize_t),
               type, url, header, ctype, body, ssl, buff, buffsize,
              )


	return nothing
end 

"""
	PetscHTTPRequest(petsclib::PetscLibType,type::Vector{Cchar}, url::Vector{Cchar}, header::Vector{Cchar}, ctype::Vector{Cchar}, body::Vector{Cchar}, sock::Cint, buff::Vector{Cchar}, buffsize::Csize_t) 
Send a request to an HTTP server

Input Parameters:
- `type`     - either "POST" or "GET"
- `url`      - URL of request host/path
- `header`   - additional header information, may be `NULL`
- `ctype`    - data type of body, for example application/json
- `body`     - data to send to server
- `sock`     - obtained with `PetscOpenSocket()`
- `buffsize` - size of buffer

Output Parameter:
- `buff` - everything returned from server

Level: advanced

-seealso: `PetscHTTPSRequest()`, `PetscOpenSocket()`, `PetscHTTPSConnect()`, `PetscPullJSONValue()`

# External Links
$(_doc_external("Sys/PetscHTTPRequest"))
"""
function PetscHTTPRequest(petsclib::PetscLibType, type::Vector{Cchar}, url::Vector{Cchar}, header::Vector{Cchar}, ctype::Vector{Cchar}, body::Vector{Cchar}, sock::Cint, buff::Vector{Cchar}, buffsize::Csize_t) end

@for_petsc function PetscHTTPRequest(petsclib::$UnionPetscLib, type::Vector{Cchar}, url::Vector{Cchar}, header::Vector{Cchar}, ctype::Vector{Cchar}, body::Vector{Cchar}, sock::Cint, buff::Vector{Cchar}, buffsize::Csize_t )

    @chk ccall(
               (:PetscHTTPRequest, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Cchar}, Csize_t),
               type, url, header, ctype, body, sock, buff, buffsize,
              )


	return nothing
end 

"""
	PetscHTTPSConnect(petsclib::PetscLibType,host::Vector{Cchar}, port::Cint, ctx::SSL_CTX, sock::Cint, ssl::SSL) 
connect to a HTTPS server

Input Parameters:
- `host` - the name of the machine hosting the HTTPS server
- `port` - the port number where the server is hosting, usually 443
- `ctx`  - value obtained with `PetscSSLInitializeContext()`

Output Parameters:
- `sock` - socket to connect
- `ssl`  - the argument passed to `PetscHTTPSRequest()`

Level: advanced

-seealso: `PetscOpenSocket()`, `PetscHTTPSRequest()`, `PetscSSLInitializeContext()`

# External Links
$(_doc_external("Sys/PetscHTTPSConnect"))
"""
function PetscHTTPSConnect(petsclib::PetscLibType, host::Vector{Cchar}, port::Cint, ctx::SSL_CTX, sock::Cint, ssl::SSL) end

@for_petsc function PetscHTTPSConnect(petsclib::$UnionPetscLib, host::Vector{Cchar}, port::Cint, ctx::SSL_CTX, sock::Cint, ssl::SSL )

    @chk ccall(
               (:PetscHTTPSConnect, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cint, Ptr{SSL_CTX}, Ptr{Cint}, SSL),
               host, port, ctx, sock, ssl,
              )


	return nothing
end 

"""
	found::PetscBool = PetscPullJSONValue(petsclib::PetscLibType,buff::Vector{Cchar}, key::Vector{Cchar}, value::Vector{Cchar}, valuelen::Csize_t) 
Given a JSON response containing the substring with "key" : "value"  where there may or not be spaces around the : returns the value.

Input Parameters:
- `buff`     - the char array containing the possible values
- `key`      - the key of the requested value
- `valuelen` - the length of the array to contain the value associated with the key

Output Parameters:
- `value` - the value obtained
- `found` - flag indicating if the value was found in the buff

Level: advanced

-seealso: `PetscOpenSocket()`, `PetscHTTPSRequest()`, `PetscSSLInitializeContext()`, `PetscPushJSONValue()`

# External Links
$(_doc_external("Sys/PetscPullJSONValue"))
"""
function PetscPullJSONValue(petsclib::PetscLibType, buff::Vector{Cchar}, key::Vector{Cchar}, value::Vector{Cchar}, valuelen::Csize_t) end

@for_petsc function PetscPullJSONValue(petsclib::$UnionPetscLib, buff::Vector{Cchar}, key::Vector{Cchar}, value::Vector{Cchar}, valuelen::Csize_t )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscPullJSONValue, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               buff, key, value, valuelen, found_,
              )

	found = found_[]

	return found
end 

"""
	PetscPushJSONValue(petsclib::PetscLibType,buff::Vector{Cchar}, key::Vector{Cchar}, value::Vector{Cchar}, bufflen::Csize_t) 
Puts a "key" : "value" pair onto a string

Input Parameters:
- `buff`    - the char array where the value will be put
- `key`     - the key value to be set
- `value`   - the value associated with the key
- `bufflen` - the size of the buffer (currently ignored)

Level: advanced

-seealso: `PetscOpenSocket()`, `PetscHTTPSRequest()`, `PetscSSLInitializeContext()`, `PetscPullJSONValue()`

# External Links
$(_doc_external("Sys/PetscPushJSONValue"))
"""
function PetscPushJSONValue(petsclib::PetscLibType, buff::Vector{Cchar}, key::Vector{Cchar}, value::Vector{Cchar}, bufflen::Csize_t) end

@for_petsc function PetscPushJSONValue(petsclib::$UnionPetscLib, buff::Vector{Cchar}, key::Vector{Cchar}, value::Vector{Cchar}, bufflen::Csize_t )

    @chk ccall(
               (:PetscPushJSONValue, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               buff, key, value, bufflen,
              )


	return nothing
end 

"""
	PetscBoxAuthorize(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, refresh_token::Vector{Cchar}, tokensize::Csize_t) 
Get authorization and refresh token for accessing Box drive from PETSc

Not Collective, only the first rank in `MPI_Comm` does anything

Input Parameters:
- `comm`      - the MPI communicator
- `tokensize` - size of the token arrays

Output Parameters:
- `access_token`  - can be used with `PetscBoxUpload()` for this one session
- `refresh_token` - can be used for ever to obtain new access_tokens with `PetscBoxRefresh()`,
guard this like a password  it gives access to your Box Drive

Level: intermediate

-seealso: `PetscBoxRefresh()`, `PetscBoxUpload()`

# External Links
$(_doc_external("Sys/PetscBoxAuthorize"))
"""
function PetscBoxAuthorize(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, refresh_token::Vector{Cchar}, tokensize::Csize_t) end

@for_petsc function PetscBoxAuthorize(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, refresh_token::Vector{Cchar}, tokensize::Csize_t )

    @chk ccall(
               (:PetscBoxAuthorize, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               comm, access_token, refresh_token, tokensize,
              )


	return nothing
end 

"""
	PetscBoxRefresh(petsclib::PetscLibType,comm::MPI_Comm, refresh_token::Vector{Cchar}, access_token::Vector{Cchar}, new_refresh_token::Vector{Cchar}, tokensize::Csize_t) 
Get a new authorization token for accessing Box drive from PETSc from a refresh token

Not Collective, only the first process in the `MPI_Comm` does anything

Input Parameters:
- `comm`          - MPI communicator
- `refresh_token` - obtained with `PetscBoxAuthorize()`, if `NULL` PETSc will first look for one in the options data
if not found it will call `PetscBoxAuthorize()`
- `tokensize`     - size of the output string access_token

Output Parameters:
- `access_token`      - token that can be passed to `PetscBoxUpload()`
- `new_refresh_token` - the old refresh token is no longer valid, not this is different than Google where the same refresh_token is used forever

Level: intermediate

-seealso: `PetscBoxAuthorize()`, `PetscBoxUpload()`

# External Links
$(_doc_external("Sys/PetscBoxRefresh"))
"""
function PetscBoxRefresh(petsclib::PetscLibType, comm::MPI_Comm, refresh_token::Vector{Cchar}, access_token::Vector{Cchar}, new_refresh_token::Vector{Cchar}, tokensize::Csize_t) end

@for_petsc function PetscBoxRefresh(petsclib::$UnionPetscLib, comm::MPI_Comm, refresh_token::Vector{Cchar}, access_token::Vector{Cchar}, new_refresh_token::Vector{Cchar}, tokensize::Csize_t )

    @chk ccall(
               (:PetscBoxRefresh, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               comm, refresh_token, access_token, new_refresh_token, tokensize,
              )


	return nothing
end 

"""
	PetscBoxUpload(petsclib::PetscLibType,comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar}) 
Loads a file to the Box Drive

This routine has not yet been written; it is just copied from Google Drive

Not collective, only the first process in the `MPI_Comm` uploads the file

Input Parameters:
- `comm`         - MPI communicator
- `access_token` - obtained with `PetscBoxRefresh()`, pass `NULL` to have PETSc generate one
- `filename`     - file to upload; if you upload multiple times it will have different names each time on Box Drive

Options Database Key:
- `-box_refresh_token XXX` - the token value

-seealso: `PetscBoxAuthorize()`, `PetscBoxRefresh()`

# External Links
$(_doc_external("Sys/PetscBoxUpload"))
"""
function PetscBoxUpload(petsclib::PetscLibType, comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar}) end

@for_petsc function PetscBoxUpload(petsclib::$UnionPetscLib, comm::MPI_Comm, access_token::Vector{Cchar}, filename::Vector{Cchar} )

    @chk ccall(
               (:PetscBoxUpload, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}),
               comm, access_token, filename,
              )


	return nothing
end 

"""
	PetscSAWsBlock(petsclib::PetscLibType) 
Blocks on SAWs until a client (person using the web browser) unblocks it

Not Collective

Level: advanced

-seealso: `PetscObjectSetName()`, `PetscObjectSAWsViewOff()`, `PetscObjectSAWsSetBlock()`, `PetscObjectSAWsBlock()`

# External Links
$(_doc_external("Sys/PetscSAWsBlock"))
"""
function PetscSAWsBlock(petsclib::PetscLibType) end

@for_petsc function PetscSAWsBlock(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSAWsBlock, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscMemoryGetCurrentUsage(petsclib::PetscLibType,mem::PetscLogDouble) 
Returns the current resident set size (memory used)
for the program.

Not Collective

Output Parameter:
- `mem` - memory usage in bytes

Options Database Key:
- `-memory_view`     - Print memory usage at end of run
- `-log_view_memory` - Display memory information for each logged event
- `-malloc_view`     - Print usage of `PetscMalloc()` in `PetscFinalize()`

Level: intermediate

-seealso: `PetscMallocGetMaximumUsage()`, `PetscMemoryGetMaximumUsage()`, `PetscMallocGetCurrentUsage()`, `PetscMemorySetGetMaximumUsage()`, `PetscMemoryView()`

# External Links
$(_doc_external("Sys/PetscMemoryGetCurrentUsage"))
"""
function PetscMemoryGetCurrentUsage(petsclib::PetscLibType, mem::PetscLogDouble) end

@for_petsc function PetscMemoryGetCurrentUsage(petsclib::$UnionPetscLib, mem::PetscLogDouble )

    @chk ccall(
               (:PetscMemoryGetCurrentUsage, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               mem,
              )


	return nothing
end 

"""
	PetscMemoryGetMaximumUsage(petsclib::PetscLibType,mem::PetscLogDouble) 
Returns the maximum resident set size (memory used)
for the program since it started (the high water mark).

Not Collective

Output Parameter:
- `mem` - memory usage in bytes

Options Database Key:
- `-memory_view`     - Print memory usage at end of run
- `-log_view_memory` - Print memory information per event
- `-malloc_view`     - Print usage of `PetscMalloc()` in `PetscFinalize()`

Level: intermediate

-seealso: `PetscMallocGetMaximumUsage()`, `PetscMemoryGetCurrentUsage()`, `PetscMallocGetCurrentUsage()`,
`PetscMemorySetGetMaximumUsage()`

# External Links
$(_doc_external("Sys/PetscMemoryGetMaximumUsage"))
"""
function PetscMemoryGetMaximumUsage(petsclib::PetscLibType, mem::PetscLogDouble) end

@for_petsc function PetscMemoryGetMaximumUsage(petsclib::$UnionPetscLib, mem::PetscLogDouble )

    @chk ccall(
               (:PetscMemoryGetMaximumUsage, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               mem,
              )


	return nothing
end 

"""
	PetscMemorySetGetMaximumUsage(petsclib::PetscLibType) 
Tells PETSc to monitor the maximum memory usage so that
`PetscMemoryGetMaximumUsage()` will work.

Not Collective

Options Database Key:
- `-memory_view`     - Print memory usage at end of run
- `-log_view_memory` - Print memory information per event
- `-malloc_view`     - Print usage of `PetscMalloc()` in `PetscFinalize()`

Level: intermediate

-seealso: `PetscMallocGetMaximumUsage()`, `PetscMemoryGetCurrentUsage()`, `PetscMallocGetCurrentUsage()`,
`PetscMemoryGetMaximumUsage()`

# External Links
$(_doc_external("Sys/PetscMemorySetGetMaximumUsage"))
"""
function PetscMemorySetGetMaximumUsage(petsclib::PetscLibType) end

@for_petsc function PetscMemorySetGetMaximumUsage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscMemorySetGetMaximumUsage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscMallocValidate(petsclib::PetscLibType,line::Cint, fnc::Vector{Cchar}, file::Vector{Cchar}) 
Test the memory for corruption.  This can be called at any time between `PetscInitialize()` and `PetscFinalize()`

Input Parameters:
- `line`     - line number where call originated.
- `function` - name of function calling
- `file`     - file where function is

Options Database Keys:
- `-malloc_test`  - turns this feature on when PETSc was not configured with `--with-debugging=0`
- `-malloc_debug` - turns this feature on anytime

Level: advanced

-seealso: `CHKMEMQ`, `PetscMalloc()`, `PetscFree()`, `PetscMallocSetDebug()`

# External Links
$(_doc_external("Sys/PetscMallocValidate"))
"""
function PetscMallocValidate(petsclib::PetscLibType, line::Cint, fnc::Vector{Cchar}, file::Vector{Cchar}) end

@for_petsc function PetscMallocValidate(petsclib::$UnionPetscLib, line::Cint, fnc::Vector{Cchar}, file::Vector{Cchar} )

    @chk ccall(
               (:PetscMallocValidate, $petsc_library),
               PetscErrorCode,
               (Cint, Ptr{Cchar}, Ptr{Cchar}),
               line, fnc, file,
              )


	return nothing
end 

"""
	PetscMemoryView(petsclib::PetscLibType,viewer::PetscViewer, message::Vector{Cchar}) 
Shows the amount of memory currently being used in a communicator.

Collective

Input Parameters:
- `viewer`  - the viewer to output the information on
- `message` - string printed before values

Options Database Keys:
- `-malloc_debug`    - have PETSc track how much memory it has allocated
- `-log_view_memory` - print memory usage per event when `-log_view` is used
- `-memory_view`     - during `PetscFinalize()` have this routine called

Level: intermediate

-seealso: `PetscMallocDump()`, `PetscMemoryGetCurrentUsage()`, `PetscMemorySetGetMaximumUsage()`, `PetscMallocView()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMemoryView"))
"""
function PetscMemoryView(petsclib::PetscLibType, viewer::PetscViewer, message::Vector{Cchar}) end

@for_petsc function PetscMemoryView(petsclib::$UnionPetscLib, viewer::PetscViewer, message::Vector{Cchar} )

    @chk ccall(
               (:PetscMemoryView, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, message,
              )


	return nothing
end 

"""
	PetscMallocGetCurrentUsage(petsclib::PetscLibType,space::PetscLogDouble) 
gets the current amount of memory used that was allocated with `PetscMalloc()`

Not Collective

Output Parameter:
- `space` - number of bytes currently allocated

Level: intermediate

-seealso: `PetscMallocDump()`, `PetscMallocGetMaximumUsage()`, `PetscMemoryGetCurrentUsage()`, `PetscMalloc()`, `PetscFree()`,
`PetscMemoryGetMaximumUsage()`

# External Links
$(_doc_external("Sys/PetscMallocGetCurrentUsage"))
"""
function PetscMallocGetCurrentUsage(petsclib::PetscLibType, space::PetscLogDouble) end

@for_petsc function PetscMallocGetCurrentUsage(petsclib::$UnionPetscLib, space::PetscLogDouble )

    @chk ccall(
               (:PetscMallocGetCurrentUsage, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               space,
              )


	return nothing
end 

"""
	PetscMallocGetMaximumUsage(petsclib::PetscLibType,space::PetscLogDouble) 
gets the maximum amount of memory used that was obtained with `PetscMalloc()` at any time
during this run, the high water mark.

Not Collective

Output Parameter:
- `space` - maximum number of bytes ever allocated at one time

Level: intermediate

-seealso: `PetscMallocDump()`, `PetscMallocView()`, `PetscMemoryGetCurrentUsage()`, `PetscMalloc()`, `PetscFree()`,
`PetscMallocPushMaximumUsage()`

# External Links
$(_doc_external("Sys/PetscMallocGetMaximumUsage"))
"""
function PetscMallocGetMaximumUsage(petsclib::PetscLibType, space::PetscLogDouble) end

@for_petsc function PetscMallocGetMaximumUsage(petsclib::$UnionPetscLib, space::PetscLogDouble )

    @chk ccall(
               (:PetscMallocGetMaximumUsage, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               space,
              )


	return nothing
end 

"""
	PetscMallocPushMaximumUsage(petsclib::PetscLibType,event::Cint) 
Adds another event to collect the maximum memory usage over an event

Not Collective

Input Parameter:
- `event` - an event id; this is just for error checking

Level: developer

-seealso: `PetscMallocDump()`, `PetscMallocView()`, `PetscMallocGetMaximumUsage()`, `PetscMemoryGetCurrentUsage()`, `PetscMalloc()`, `PetscFree()`,
`PetscMallocPopMaximumUsage()`

# External Links
$(_doc_external("Sys/PetscMallocPushMaximumUsage"))
"""
function PetscMallocPushMaximumUsage(petsclib::PetscLibType, event::Cint) end

@for_petsc function PetscMallocPushMaximumUsage(petsclib::$UnionPetscLib, event::Cint )

    @chk ccall(
               (:PetscMallocPushMaximumUsage, $petsc_library),
               PetscErrorCode,
               (Cint,),
               event,
              )


	return nothing
end 

"""
	PetscMallocPopMaximumUsage(petsclib::PetscLibType,event::Cint, mu::PetscLogDouble) 
collect the maximum memory usage over an event

Not Collective

Input Parameter:
- `event` - an event id; this is just for error checking

Output Parameter:
- `mu` - maximum amount of memory malloced during this event; high water mark relative to the beginning of the event

Level: developer

-seealso: `PetscMallocDump()`, `PetscMallocView()`, `PetscMallocGetMaximumUsage()`, `PetscMemoryGetCurrentUsage()`, `PetscMalloc()`, `PetscFree()`,
`PetscMallocPushMaximumUsage()`

# External Links
$(_doc_external("Sys/PetscMallocPopMaximumUsage"))
"""
function PetscMallocPopMaximumUsage(petsclib::PetscLibType, event::Cint, mu::PetscLogDouble) end

@for_petsc function PetscMallocPopMaximumUsage(petsclib::$UnionPetscLib, event::Cint, mu::PetscLogDouble )

    @chk ccall(
               (:PetscMallocPopMaximumUsage, $petsc_library),
               PetscErrorCode,
               (Cint, Ptr{PetscLogDouble}),
               event, mu,
              )


	return nothing
end 

"""
	PetscMallocDump(petsclib::PetscLibType,fp::Libc.FILE) 
Dumps the currently allocated memory blocks to a file. The information
printed is: size of space (in bytes), address of space, id of space,
file in which space was allocated, and line number at which it was
allocated.

Not Collective

Input Parameter:
- `fp` - file pointer.  If `fp` is `NULL`, `stdout` is assumed.

Options Database Key:
- `-malloc_dump <optional filename>` - Print summary of unfreed memory during call to `PetscFinalize()`, writing to filename if given

Level: intermediate

-seealso: `PetscMallocGetCurrentUsage()`, `PetscMallocView()`, `PetscMallocViewSet()`, `PetscMallocValidate()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocDump"))
"""
function PetscMallocDump(petsclib::PetscLibType, fp::Libc.FILE) end

@for_petsc function PetscMallocDump(petsclib::$UnionPetscLib, fp::Libc.FILE )

    @chk ccall(
               (:PetscMallocDump, $petsc_library),
               PetscErrorCode,
               (Ptr{Libc.FILE},),
               fp,
              )


	return nothing
end 

"""
	PetscMallocViewSet(petsclib::PetscLibType,logmin::PetscLogDouble) 
Activates logging of all calls to `PetscMalloc()` with a minimum size to view

Not Collective

Input Parameter:
- `logmin` - minimum allocation size to log, or `PETSC_DEFAULT` to log all memory allocations

Options Database Keys:
- `-malloc_view <optional filename>` - Activates `PetscMallocView()` in `PetscFinalize()`
- `-malloc_view_threshold <min>`     - Sets a minimum size if `-malloc_view` is used
- `-log_view_memory`                 - view the memory usage also with the -log_view option

Level: advanced

-seealso: `PetscMallocViewGet()`, `PetscMallocDump()`, `PetscMallocView()`, `PetscMallocTraceSet()`, `PetscMallocValidate()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocViewSet"))
"""
function PetscMallocViewSet(petsclib::PetscLibType, logmin::PetscLogDouble) end

@for_petsc function PetscMallocViewSet(petsclib::$UnionPetscLib, logmin::PetscLogDouble )

    @chk ccall(
               (:PetscMallocViewSet, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               logmin,
              )


	return nothing
end 

"""
	logging::PetscBool = PetscMallocViewGet(petsclib::PetscLibType) 
Determine whether calls to `PetscMalloc()` are being logged

Not Collective

Output Parameter:
- `logging` - `PETSC_TRUE` if logging is active

Options Database Key:
- `-malloc_view <optional filename>` - Activates `PetscMallocView()`

Level: advanced

-seealso: `PetscMallocViewSet()`, `PetscMallocDump()`, `PetscMallocView()`, `PetscMallocTraceGet()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocViewGet"))
"""
function PetscMallocViewGet(petsclib::PetscLibType) end

@for_petsc function PetscMallocViewGet(petsclib::$UnionPetscLib)
	logging_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscMallocViewGet, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               logging_,
              )

	logging = logging_[]

	return logging
end 

"""
	PetscMallocTraceSet(petsclib::PetscLibType,viewer::PetscViewer, active::PetscBool, logmin::PetscLogDouble) 
Trace all calls to `PetscMalloc()`. That is print each `PetscMalloc()` and `PetscFree()` call to a viewer.

Not Collective

Input Parameters:
- `viewer` - The viewer to use for tracing, or `NULL` to use `PETSC_VIEWER_STDOUT_SELF`
- `active` - Flag to activate or deactivate tracing
- `logmin` - The smallest memory size that will be logged

Level: advanced

-seealso: `PetscMallocTraceGet()`, `PetscMallocViewGet()`, `PetscMallocDump()`, `PetscMallocView()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocTraceSet"))
"""
function PetscMallocTraceSet(petsclib::PetscLibType, viewer::PetscViewer, active::PetscBool, logmin::PetscLogDouble) end

@for_petsc function PetscMallocTraceSet(petsclib::$UnionPetscLib, viewer::PetscViewer, active::PetscBool, logmin::PetscLogDouble )

    @chk ccall(
               (:PetscMallocTraceSet, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool, PetscLogDouble),
               viewer, active, logmin,
              )


	return nothing
end 

"""
	logging::PetscBool = PetscMallocTraceGet(petsclib::PetscLibType) 
Determine whether all calls to `PetscMalloc()` are being traced

Not Collective

Output Parameter:
- `logging` - `PETSC_TRUE` if logging is active

Options Database Key:
- `-malloc_view <optional filename>` - Activates `PetscMallocView()`

Level: advanced

This only does anything if `-malloc_debug` (or `-malloc_test` if PETSc was configured with debugging) has been used

-seealso: `PetscMallocTraceSet()`, `PetscMallocViewGet()`, `PetscMallocDump()`, `PetscMallocView()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocTraceGet"))
"""
function PetscMallocTraceGet(petsclib::PetscLibType) end

@for_petsc function PetscMallocTraceGet(petsclib::$UnionPetscLib)
	logging_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscMallocTraceGet, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               logging_,
              )

	logging = logging_[]

	return logging
end 

"""
	PetscMallocView(petsclib::PetscLibType,fp::Libc.FILE) 
Saves the log of all calls to `PetscMalloc()`; also calls `PetscMemoryGetMaximumUsage()`

Not Collective

Input Parameter:
- `fp` - file pointer; or `NULL`

Options Database Key:
- `-malloc_view <optional filename>` - Activates `PetscMallocView()` in `PetscFinalize()`

Level: advanced

-seealso: `PetscMallocGetCurrentUsage()`, `PetscMallocDump()`, `PetscMallocViewSet()`, `PetscMemoryView()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocView"))
"""
function PetscMallocView(petsclib::PetscLibType, fp::Libc.FILE) end

@for_petsc function PetscMallocView(petsclib::$UnionPetscLib, fp::Libc.FILE )

    @chk ccall(
               (:PetscMallocView, $petsc_library),
               PetscErrorCode,
               (Ptr{Libc.FILE},),
               fp,
              )


	return nothing
end 

"""
	PetscMallocSetDebug(petsclib::PetscLibType,eachcall::PetscBool, initializenan::PetscBool) 
Set's PETSc memory debugging

Not Collective

Input Parameters:
- `eachcall`      - checks the entire heap of allocated memory for issues on each call to `PetscMalloc()` and `PetscFree()`, slow
- `initializenan` - initializes all memory with `NaN` to catch use of uninitialized floating point arrays

Options Database Keys:
- `-malloc_debug <true or false>` - turns on or off debugging
- `-malloc_test`                  - turns on all debugging if PETSc was configured with debugging including `-malloc_dump`, otherwise ignored
- `-malloc_view_threshold t`      - log only allocations larger than t
- `-malloc_dump <filename>`       - print a list of all memory that has not been freed, in `PetscFinalize()`

Level: developer

-seealso: `CHKMEMQ`, `PetscMallocValidate()`, `PetscMallocGetDebug()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocSetDebug"))
"""
function PetscMallocSetDebug(petsclib::PetscLibType, eachcall::PetscBool, initializenan::PetscBool) end

@for_petsc function PetscMallocSetDebug(petsclib::$UnionPetscLib, eachcall::PetscBool, initializenan::PetscBool )

    @chk ccall(
               (:PetscMallocSetDebug, $petsc_library),
               PetscErrorCode,
               (PetscBool, PetscBool),
               eachcall, initializenan,
              )


	return nothing
end 

"""
	basic::PetscBool,eachcall::PetscBool,initializenan::PetscBool = PetscMallocGetDebug(petsclib::PetscLibType) 
Indicates what PETSc memory debugging it is doing.

Not Collective

Output Parameters:
- `basic`         - doing basic debugging
- `eachcall`      - checks the entire memory heap at each `PetscMalloc()`/`PetscFree()`
- `initializenan` - initializes memory with `NaN`

Level: intermediate

-seealso: `CHKMEMQ`, `PetscMallocValidate()`, `PetscMallocSetDebug()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocGetDebug"))
"""
function PetscMallocGetDebug(petsclib::PetscLibType) end

@for_petsc function PetscMallocGetDebug(petsclib::$UnionPetscLib)
	basic_ = Ref{PetscBool}()
	eachcall_ = Ref{PetscBool}()
	initializenan_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscMallocGetDebug, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}),
               basic_, eachcall_, initializenan_,
              )

	basic = basic_[]
	eachcall = eachcall_[]
	initializenan = initializenan_[]

	return basic,eachcall,initializenan
end 

"""
	PetscMallocLogRequestedSizeSet(petsclib::PetscLibType,flg::PetscBool) 
Whether to log the requested or aligned memory size

Not Collective

Input Parameter:
- `flg` - `PETSC_TRUE` to log the requested memory size

Options Database Key:
- `-malloc_requested_size <bool>` - Sets this flag

Level: developer

-seealso: `PetscMallocLogRequestedSizeGet()`, `PetscMallocViewSet()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocLogRequestedSizeSet"))
"""
function PetscMallocLogRequestedSizeSet(petsclib::PetscLibType, flg::PetscBool) end

@for_petsc function PetscMallocLogRequestedSizeSet(petsclib::$UnionPetscLib, flg::PetscBool )

    @chk ccall(
               (:PetscMallocLogRequestedSizeSet, $petsc_library),
               PetscErrorCode,
               (PetscBool,),
               flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscMallocLogRequestedSizeGet(petsclib::PetscLibType) 
Whether to log the requested or aligned memory size

Not Collective

Output Parameter:
- `flg` - `PETSC_TRUE` if we log the requested memory size

Level: developer

-seealso: `PetscMallocLogRequestedSizeSet()`, `PetscMallocViewSet()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocLogRequestedSizeGet"))
"""
function PetscMallocLogRequestedSizeGet(petsclib::PetscLibType) end

@for_petsc function PetscMallocLogRequestedSizeGet(petsclib::$UnionPetscLib)
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscMallocLogRequestedSizeGet, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscMallocSet(petsclib::PetscLibType,imalloc::external, ifree::external, iralloc::external) 
Sets the underlying allocation routines used by `PetscMalloc()` and `PetscFree()`

Not Collective, No Fortran Support

Input Parameters:
- `imalloc` - the routine that provides the `malloc()` implementation (also provides `calloc()`, which is used depending on the second argument)
- `ifree`   - the routine that provides the `free()` implementation
- `iralloc` - the routine that provides the `realloc()` implementation

Level: developer

-seealso: `PetscMallocClear()`, `PetscInitialize()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocSet"))
"""
function PetscMallocSet(petsclib::PetscLibType, imalloc::external, ifree::external, iralloc::external) end

@for_petsc function PetscMallocSet(petsclib::$UnionPetscLib, imalloc::external, ifree::external, iralloc::external )

    @chk ccall(
               (:PetscMallocSet, $petsc_library),
               PetscErrorCode,
               (external, external, external),
               imalloc, ifree, iralloc,
              )


	return nothing
end 

"""
	PetscMallocClear(petsclib::PetscLibType) 
Resets the routines used by `PetscMalloc()` and `PetscFree()`

Not Collective

Level: developer

-seealso: `PetscMallocSet()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocClear"))
"""
function PetscMallocClear(petsclib::PetscLibType) end

@for_petsc function PetscMallocClear(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscMallocClear, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscMemoryTrace(petsclib::PetscLibType,label::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscMemoryTrace"))
"""
function PetscMemoryTrace(petsclib::PetscLibType, label::Vector{Cchar}) end

@for_petsc function PetscMemoryTrace(petsclib::$UnionPetscLib, label::Vector{Cchar} )

    @chk ccall(
               (:PetscMemoryTrace, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               label,
              )


	return nothing
end 

"""
	PetscMallocSetDRAM(petsclib::PetscLibType) 
Set `PetscMalloc()` to use DRAM.
If memkind is available, change the memkind type. Otherwise, switch the
current malloc and free routines to the `PetscMallocAlign()` and
`PetscFreeAlign()` (PETSc default).

Not Collective

Level: developer

-seealso: `PetscMallocReset()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocSetDRAM"))
"""
function PetscMallocSetDRAM(petsclib::PetscLibType) end

@for_petsc function PetscMallocSetDRAM(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscMallocSetDRAM, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscMallocResetDRAM(petsclib::PetscLibType) 
Reset the changes made by `PetscMallocSetDRAM()`

Not Collective

Level: developer

-seealso: `PetscMallocSetDRAM()`

# External Links
$(_doc_external("Sys/PetscMallocResetDRAM"))
"""
function PetscMallocResetDRAM(petsclib::PetscLibType) end

@for_petsc function PetscMallocResetDRAM(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscMallocResetDRAM, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscMallocSetCoalesce(petsclib::PetscLibType,coalesce::PetscBool) 
Use coalesced `PetscMalloc()` when allocating groups of objects, that is when using `PetscMallocN()`

Not Collective

Input Parameter:
- `coalesce` - `PETSC_TRUE` to use coalesced malloc for multi-memory allocation.

Options Database Key:
- `-malloc_coalesce` - turn coalesced `PetscMallocN()` on or off

Level: developer

-seealso: `PetscMallocA()`, `PetscMalloc()`, `PetscFree()`

# External Links
$(_doc_external("Sys/PetscMallocSetCoalesce"))
"""
function PetscMallocSetCoalesce(petsclib::PetscLibType, coalesce::PetscBool) end

@for_petsc function PetscMallocSetCoalesce(petsclib::$UnionPetscLib, coalesce::PetscBool )

    @chk ccall(
               (:PetscMallocSetCoalesce, $petsc_library),
               PetscErrorCode,
               (PetscBool,),
               coalesce,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscTestFile(petsclib::PetscLibType,fname::Vector{Cchar}, mode::Cchar) 
checks for the existence of a file

Not Collective

Input Parameters:
- `fname` - the filename
- `mode`  - either 'r', 'w', 'x' or '\0'

Output Parameter:
- `flg` - the file exists and satisfies the mode

Level: intermediate

-seealso: `PetscTestDirectory()`, `PetscLs()`

# External Links
$(_doc_external("Sys/PetscTestFile"))
"""
function PetscTestFile(petsclib::PetscLibType, fname::Vector{Cchar}, mode::Cchar) end

@for_petsc function PetscTestFile(petsclib::$UnionPetscLib, fname::Vector{Cchar}, mode::Cchar )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscTestFile, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{PetscBool}),
               fname, mode, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PetscTestDirectory(petsclib::PetscLibType,dirname::Vector{Cchar}, mode::Cchar) 
checks for the existence of a directory

Not Collective

Input Parameters:
- `dirname` - the directory name
- `mode`    - either 'r', 'w', or 'x'

Output Parameter:
- `flg` - the directory exists and satisfies the mode

Level: intermediate

-seealso: `PetscTestFile()`, `PetscLs()`, `PetscRMTree()`

# External Links
$(_doc_external("Sys/PetscTestDirectory"))
"""
function PetscTestDirectory(petsclib::PetscLibType, dirname::Vector{Cchar}, mode::Cchar) end

@for_petsc function PetscTestDirectory(petsclib::$UnionPetscLib, dirname::Vector{Cchar}, mode::Cchar )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscTestDirectory, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{PetscBool}),
               dirname, mode, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PetscLs(petsclib::PetscLibType,comm::MPI_Comm, dirname::Vector{Cchar}, found::Vector{Cchar}, tlen::Csize_t) 
produce a listing of the files in a directory

Collective

Input Parameters:
- `comm`    - the MPI communicator
- `dirname` - the directory name
- `tlen`    - the length of the buffer `found`

Output Parameters:
- `found` - listing of files
- `flg`   - the directory exists

Level: intermediate

-seealso: `PetscTestFile()`, `PetscRMTree()`, `PetscTestDirectory()`

# External Links
$(_doc_external("Sys/PetscLs"))
"""
function PetscLs(petsclib::PetscLibType, comm::MPI_Comm, dirname::Vector{Cchar}, found::Vector{Cchar}, tlen::Csize_t) end

@for_petsc function PetscLs(petsclib::$UnionPetscLib, comm::MPI_Comm, dirname::Vector{Cchar}, found::Vector{Cchar}, tlen::Csize_t )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLs, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               comm, dirname, found, tlen, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscGetRealPath(petsclib::PetscLibType,path::Vector{Cchar}, rpath::Vector{Cchar}) 
Get the path without symbolic links etc. in absolute form.

Not Collective

Input Parameter:
- `path` - path to resolve

Output Parameter:
- `rpath` - resolved path

Level: developer

-seealso: `PetscGetFullPath()`

# External Links
$(_doc_external("Sys/PetscGetRealPath"))
"""
function PetscGetRealPath(petsclib::PetscLibType, path::Vector{Cchar}, rpath::Vector{Cchar}) end

@for_petsc function PetscGetRealPath(petsclib::$UnionPetscLib, path::Vector{Cchar}, rpath::Vector{Cchar} )

    @chk ccall(
               (:PetscGetRealPath, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               path, rpath,
              )


	return nothing
end 

"""
	PetscGetTmp(petsclib::PetscLibType,comm::MPI_Comm, dir::Vector{Cchar}, len::Csize_t) 
Gets the name of the "tmp" directory, often this is `/tmp`

Collective

Input Parameters:
- `comm` - MPI_Communicator that may share tmp
- `len`  - length of string to hold name

Output Parameter:
- `dir` - directory name

Options Database Keys:
- `-shared_tmp`     - indicates the directory is known to be shared among the MPI processes
- `-not_shared_tmp` - indicates the directory is known to be not shared among the MPI processes
- `-tmp tmpdir`     - name of the directory you wish to use as tmp

Environmental Variables:
- `PETSC_SHARED_TMP`     - indicates the directory is known to be shared among the MPI processes
- `PETSC_NOT_SHARED_TMP` - indicates the directory is known to be not shared among the MPI processes
- `PETSC_TMP`            - name of the directory you wish to use as tmp

Level: developer

-seealso: `PetscSharedTmp()`, `PetscSharedWorkingDirectory()`, `PetscGetWorkingDirectory()`, `PetscGetHomeDirectory()`

# External Links
$(_doc_external("Sys/PetscGetTmp"))
"""
function PetscGetTmp(petsclib::PetscLibType, comm::MPI_Comm, dir::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscGetTmp(petsclib::$UnionPetscLib, comm::MPI_Comm, dir::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscGetTmp, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Csize_t),
               comm, dir, len,
              )


	return nothing
end 

"""
	shared::PetscBool = PetscSharedTmp(petsclib::PetscLibType,comm::MPI_Comm) 
Determines if all processors in a communicator share a
tmp directory or have different ones.

Collective

Input Parameter:
- `comm` - MPI_Communicator that may share tmp

Output Parameter:
- `shared` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
- `-shared_tmp`     - indicates the directory is known to be shared among the MPI processes
- `-not_shared_tmp` - indicates the directory is known to be not shared among the MPI processes
- `-tmp tmpdir`     - name of the directory you wish to use as tmp

Environmental Variables:
- `PETSC_SHARED_TMP`     - indicates the directory is known to be shared among the MPI processes
- `PETSC_NOT_SHARED_TMP` - indicates the directory is known to be not shared among the MPI processes
- `PETSC_TMP`            - name of the directory you wish to use as tmp

Level: developer

-seealso: `PetscGetTmp()`, `PetscSharedWorkingDirectory()`, `PetscGetWorkingDirectory()`, `PetscGetHomeDirectory()`

# External Links
$(_doc_external("Sys/PetscSharedTmp"))
"""
function PetscSharedTmp(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscSharedTmp(petsclib::$UnionPetscLib, comm::MPI_Comm )
	shared_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSharedTmp, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscBool}),
               comm, shared_,
              )

	shared = shared_[]

	return shared
end 

"""
	shared::PetscBool = PetscSharedWorkingDirectory(petsclib::PetscLibType,comm::MPI_Comm) 
Determines if all processors in a communicator share a working directory or have different ones.

Collective

Input Parameter:
- `comm` - MPI_Communicator that may share working directory

Output Parameter:
- `shared` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
- `-shared_working_directory`     - indicates the directory is known to be shared among the MPI processes
- `-not_shared_working_directory` - indicates the directory is known to be not shared among the MPI processes

Environmental Variables:
- `PETSC_SHARED_WORKING_DIRECTORY`     - indicates the directory is known to be shared among the MPI processes
- `PETSC_NOT_SHARED_WORKING_DIRECTORY` - indicates the directory is known to be not shared among the MPI processes

Level: developer

-seealso: `PetscGetTmp()`, `PetscSharedTmp()`, `PetscGetWorkingDirectory()`, `PetscGetHomeDirectory()`

# External Links
$(_doc_external("Sys/PetscSharedWorkingDirectory"))
"""
function PetscSharedWorkingDirectory(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscSharedWorkingDirectory(petsclib::$UnionPetscLib, comm::MPI_Comm )
	shared_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSharedWorkingDirectory, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscBool}),
               comm, shared_,
              )

	shared = shared_[]

	return shared
end 

"""
	found::PetscBool = PetscFileRetrieve(petsclib::PetscLibType,comm::MPI_Comm, url::Vector{Cchar}, locname::Vector{Cchar}, llen::Csize_t) 
Obtains a file from a URL or a compressed file
and copies into local disk space as uncompressed.

Collective

Input Parameters:
- `comm` - processors accessing the file
- `url`  - name of file, including entire URL (with or without .gz)
- `llen` - length of `localname`

Output Parameters:
- `localname` - name of local copy of file - valid on only process zero
- `found`     - if found or retrieved the file - valid on all processes

Level: developer

-seealso: `PetscDLLibraryRetrieve()`

# External Links
$(_doc_external("Sys/PetscFileRetrieve"))
"""
function PetscFileRetrieve(petsclib::PetscLibType, comm::MPI_Comm, url::Vector{Cchar}, locname::Vector{Cchar}, llen::Csize_t) end

@for_petsc function PetscFileRetrieve(petsclib::$UnionPetscLib, comm::MPI_Comm, url::Vector{Cchar}, locname::Vector{Cchar}, llen::Csize_t )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscFileRetrieve, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               comm, url, locname, llen, found_,
              )

	found = found_[]

	return found
end 

"""
	PetscGetFullPath(petsclib::PetscLibType,path::Vector{Cchar}, fullpath::Vector{Cchar}, flen::Csize_t) 
Given a filename, returns the fully qualified file name.

Not Collective

Input Parameters:
- `path` - pathname to qualify
- `flen` - size of `fullpath`

Output Parameter:
- `fullpath` - buffer to hold the full pathname

Level: developer

-seealso: `PetscGetRelativePath()`

# External Links
$(_doc_external("Sys/PetscGetFullPath"))
"""
function PetscGetFullPath(petsclib::PetscLibType, path::Vector{Cchar}, fullpath::Vector{Cchar}, flen::Csize_t) end

@for_petsc function PetscGetFullPath(petsclib::$UnionPetscLib, path::Vector{Cchar}, fullpath::Vector{Cchar}, flen::Csize_t )

    @chk ccall(
               (:PetscGetFullPath, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               path, fullpath, flen,
              )


	return nothing
end 

"""
	PetscByteSwap(petsclib::PetscLibType,data::Cvoid, pdtype::PetscDataType, count::PetscCount) 

# External Links
$(_doc_external("Sys/PetscByteSwap"))
"""
function PetscByteSwap(petsclib::PetscLibType, data::Cvoid, pdtype::PetscDataType, count::PetscCount) end

@for_petsc function PetscByteSwap(petsclib::$UnionPetscLib, data::Cvoid, pdtype::PetscDataType, count::PetscCount )

    @chk ccall(
               (:PetscByteSwap, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, PetscDataType, PetscCount),
               data, pdtype, count,
              )


	return nothing
end 

"""
	count::PetscInt = PetscBinaryRead(petsclib::PetscLibType,fd::Cint, data::Cvoid, num::PetscCount, type::PetscDataType) 
Reads from a binary file.

Not Collective

Input Parameters:
- `fd`   - the file descriptor
- `num`  - the maximum number of items to read
- `type` - the type of items to read (`PETSC_INT`, `PETSC_REAL`, `PETSC_SCALAR`, etc.)

Output Parameters:
- `data`  - the buffer, this is an array of the type that matches the value in `type`
- `count` - the number of items read, optional

Level: developer

-seealso: `PetscBinaryWrite()`, `PetscBinaryOpen()`, `PetscBinaryClose()`, `PetscViewerBinaryGetDescriptor()`, `PetscBinarySynchronizedWrite()`,
`PetscBinarySynchronizedRead()`, `PetscBinarySynchronizedSeek()`

# External Links
$(_doc_external("Sys/PetscBinaryRead"))
"""
function PetscBinaryRead(petsclib::PetscLibType, fd::Cint, data::Cvoid, num::PetscCount, type::PetscDataType) end

@for_petsc function PetscBinaryRead(petsclib::$UnionPetscLib, fd::Cint, data::Cvoid, num::PetscCount, type::PetscDataType )
	count_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscBinaryRead, $petsc_library),
               PetscErrorCode,
               (Cint, Ptr{Cvoid}, PetscCount, Ptr{$PetscInt}, PetscDataType),
               fd, data, num, count_, type,
              )

	count = count_[]

	return count
end 

"""
	PetscBinaryWrite(petsclib::PetscLibType,fd::Cint, p::Cvoid, n::PetscCount, type::PetscDataType) 
Writes to a binary file.

Not Collective

Input Parameters:
- `fd`   - the file
- `p`    - the buffer, an array of the type that matches the value in `type`
- `n`    - the number of items to write
- `type` - the type of items to read (`PETSC_INT`, `PETSC_REAL` or `PETSC_SCALAR`)

Level: advanced

-seealso: `PetscBinaryRead()`, `PetscBinaryOpen()`, `PetscBinaryClose()`, `PetscViewerBinaryGetDescriptor()`, `PetscBinarySynchronizedWrite()`,
`PetscBinarySynchronizedRead()`, `PetscBinarySynchronizedSeek()`

# External Links
$(_doc_external("Sys/PetscBinaryWrite"))
"""
function PetscBinaryWrite(petsclib::PetscLibType, fd::Cint, p::Cvoid, n::PetscCount, type::PetscDataType) end

@for_petsc function PetscBinaryWrite(petsclib::$UnionPetscLib, fd::Cint, p::Cvoid, n::PetscCount, type::PetscDataType )

    @chk ccall(
               (:PetscBinaryWrite, $petsc_library),
               PetscErrorCode,
               (Cint, Ptr{Cvoid}, PetscCount, PetscDataType),
               fd, p, n, type,
              )


	return nothing
end 

"""
	PetscBinaryOpen(petsclib::PetscLibType,name::Vector{Cchar}, mode::PetscFileMode, fd::Cint) 
Opens a PETSc binary file.

Not Collective

Input Parameters:
- `name` - filename
- `mode` - open mode of binary file, one of `FILE_MODE_READ`, `FILE_MODE_WRITE`, `FILE_MODE_APPEND`

Output Parameter:
- `fd` - the file

Level: advanced

-seealso: `PetscBinaryRead()`, `PetscBinaryWrite()`, `PetscFileMode`, `PetscViewerFileSetMode()`, `PetscViewerBinaryGetDescriptor()`,
`PetscBinarySynchronizedWrite()`, `PetscBinarySynchronizedRead()`, `PetscBinarySynchronizedSeek()`

# External Links
$(_doc_external("Sys/PetscBinaryOpen"))
"""
function PetscBinaryOpen(petsclib::PetscLibType, name::Vector{Cchar}, mode::PetscFileMode, fd::Cint) end

@for_petsc function PetscBinaryOpen(petsclib::$UnionPetscLib, name::Vector{Cchar}, mode::PetscFileMode, fd::Cint )

    @chk ccall(
               (:PetscBinaryOpen, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscFileMode, Ptr{Cint}),
               name, mode, fd,
              )


	return nothing
end 

"""
	PetscBinaryClose(petsclib::PetscLibType,fd::Cint) 
Closes a PETSc binary file.

Not Collective

Output Parameter:
- `fd` - the file

Level: advanced

-seealso: `PetscBinaryRead()`, `PetscBinaryWrite()`, `PetscBinaryOpen()`, `PetscBinarySynchronizedWrite()`, `PetscBinarySynchronizedRead()`,
`PetscBinarySynchronizedSeek()`

# External Links
$(_doc_external("Sys/PetscBinaryClose"))
"""
function PetscBinaryClose(petsclib::PetscLibType, fd::Cint) end

@for_petsc function PetscBinaryClose(petsclib::$UnionPetscLib, fd::Cint )

    @chk ccall(
               (:PetscBinaryClose, $petsc_library),
               PetscErrorCode,
               (Cint,),
               fd,
              )


	return nothing
end 

"""
	count::PetscInt = PetscBinarySynchronizedRead(petsclib::PetscLibType,comm::MPI_Comm, fd::Cint, data::Cvoid, num::PetscInt, type::PetscDataType) 
Reads from a binary file, all MPI processes get the same values

Collective, No Fortran Support

Input Parameters:
- `comm` - the MPI communicator
- `fd`   - the file descriptor
- `num`  - the maximum number of items to read
- `type` - the type of items to read (`PETSC_INT`, `PETSC_REAL`, `PETSC_SCALAR`, etc.)

Output Parameters:
- `data`  - the buffer, an array of the type that matches the value in `type`
- `count` - the number of items read, optional

Level: developer

-seealso: `PetscBinaryWrite()`, `PetscBinaryOpen()`, `PetscBinaryClose()`, `PetscBinaryRead()`, `PetscBinarySynchronizedWrite()`,
`PetscBinarySynchronizedSeek()`

# External Links
$(_doc_external("Sys/PetscBinarySynchronizedRead"))
"""
function PetscBinarySynchronizedRead(petsclib::PetscLibType, comm::MPI_Comm, fd::Cint, data::Cvoid, num::PetscInt, type::PetscDataType) end

@for_petsc function PetscBinarySynchronizedRead(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Cint, data::Cvoid, num::$PetscInt, type::PetscDataType )
	count_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscBinarySynchronizedRead, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Cint, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
               comm, fd, data, num, count_, type,
              )

	count = count_[]

	return count
end 

"""
	PetscBinarySynchronizedWrite(petsclib::PetscLibType,comm::MPI_Comm, fd::Cint, p::Cvoid, n::PetscInt, type::PetscDataType) 
writes to a binary file.

Collective, No Fortran Support

Input Parameters:
- `comm` - the MPI communicator
- `fd`   - the file
- `n`    - the number of items to write
- `p`    - the buffer, an array of the type that matches the value in `type`
- `type` - the type of items to write (`PETSC_INT`, `PETSC_REAL` or `PETSC_SCALAR`)

Level: developer

-seealso: `PetscBinaryWrite()`, `PetscBinaryOpen()`, `PetscBinaryClose()`, `PetscBinaryRead()`, `PetscBinarySynchronizedRead()`,
`PetscBinarySynchronizedSeek()`

# External Links
$(_doc_external("Sys/PetscBinarySynchronizedWrite"))
"""
function PetscBinarySynchronizedWrite(petsclib::PetscLibType, comm::MPI_Comm, fd::Cint, p::Cvoid, n::PetscInt, type::PetscDataType) end

@for_petsc function PetscBinarySynchronizedWrite(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Cint, p::Cvoid, n::$PetscInt, type::PetscDataType )

    @chk ccall(
               (:PetscBinarySynchronizedWrite, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Cint, Ptr{Cvoid}, $PetscInt, PetscDataType),
               comm, fd, p, n, type,
              )


	return nothing
end 

"""
	PetscStartMatlab(petsclib::PetscLibType,comm::MPI_Comm, machine::Vector{Cchar}, script::Vector{Cchar}, fp::Libc.FILE) 
starts up MATLAB with a MATLAB script

Logically Collective, but only MPI rank 0 in the communicator does anything

Input Parameters:
- `comm`    - MPI communicator
- `machine` - optional machine to run MATLAB on
- `script`  - name of script (without the .m)

Output Parameter:
- `fp` - a file pointer returned from `PetscPOpen()`

Level: intermediate

-seealso: `PetscPOpen()`, `PetscPClose()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscStartMatlab"))
"""
function PetscStartMatlab(petsclib::PetscLibType, comm::MPI_Comm, machine::Vector{Cchar}, script::Vector{Cchar}, fp::Libc.FILE) end

@for_petsc function PetscStartMatlab(petsclib::$UnionPetscLib, comm::MPI_Comm, machine::Vector{Cchar}, script::Vector{Cchar}, fp::Libc.FILE )

    @chk ccall(
               (:PetscStartMatlab, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Libc.FILE),
               comm, machine, script, fp,
              )


	return nothing
end 

"""
	PetscFOpen(petsclib::PetscLibType,comm::MPI_Comm, name::Vector{Cchar}, mode::Vector{Cchar}, fp::Libc.FILE) 
Has the first process in the MPI communicator open a file;
all others do nothing.

Logically Collective

Input Parameters:
- `comm` - the MPI communicator
- `name` - the filename
- `mode` - the mode for `fopen()`, usually "w"

Output Parameter:
- `fp` - the file pointer

Level: developer

-seealso: `PetscFClose()`, `PetscSynchronizedFGets()`, `PetscSynchronizedPrintf()`, `PetscSynchronizedFlush()`,
`PetscFPrintf()`

# External Links
$(_doc_external("Sys/PetscFOpen"))
"""
function PetscFOpen(petsclib::PetscLibType, comm::MPI_Comm, name::Vector{Cchar}, mode::Vector{Cchar}, fp::Libc.FILE) end

@for_petsc function PetscFOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::Vector{Cchar}, mode::Vector{Cchar}, fp::Libc.FILE )

    @chk ccall(
               (:PetscFOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Libc.FILE),
               comm, name, mode, fp,
              )


	return nothing
end 

"""
	PetscFClose(petsclib::PetscLibType,comm::MPI_Comm, fd::Libc.FILE) 
Has MPI rank 0 in the communicator close a
file (usually obtained with `PetscFOpen()`; all others do nothing.

Logically Collective

Input Parameters:
- `comm` - the MPI communicator
- `fd`   - the file, opened with `PetscFOpen()`

Level: developer

-seealso: `PetscFOpen()`

# External Links
$(_doc_external("Sys/PetscFClose"))
"""
function PetscFClose(petsclib::PetscLibType, comm::MPI_Comm, fd::Libc.FILE) end

@for_petsc function PetscFClose(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Libc.FILE )

    @chk ccall(
               (:PetscFClose, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}),
               comm, fd,
              )


	return nothing
end 

"""
	PetscPClose(petsclib::PetscLibType,comm::MPI_Comm, fd::Libc.FILE) 
Closes (ends) a program on MPI rank 0 run with `PetscPOpen()`

Collective, but only MPI rank 0 does anything

Input Parameters:
- `comm` - MPI communicator, only rank 0 performs the close
- `fd`   - the file pointer where program input or output may be read or `NULL` if don't care

Level: intermediate

-seealso: `PetscFOpen()`, `PetscFClose()`, `PetscPOpen()`

# External Links
$(_doc_external("Sys/PetscPClose"))
"""
function PetscPClose(petsclib::PetscLibType, comm::MPI_Comm, fd::Libc.FILE) end

@for_petsc function PetscPClose(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Libc.FILE )

    @chk ccall(
               (:PetscPClose, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}),
               comm, fd,
              )


	return nothing
end 

"""
	PetscPOpen(petsclib::PetscLibType,comm::MPI_Comm, machine::Vector{Cchar}, program::Vector{Cchar}, mode::Vector{Cchar}, fp::Libc.FILE) 
Runs a program on MPI rank 0 and sends either its input or output to
a file.

Logically Collective, but only MPI rank 0 runs the command

Input Parameters:
- `comm`    - MPI communicator, only processor zero runs the program
- `machine` - machine to run command on or `NULL`, or a string with 0 in first location
- `program` - name of program to run
- `mode`    - either "r" or "w"

Output Parameter:
- `fp` - the file pointer where program input or output may be read or `NULL` if results are not needed

Level: intermediate

-seealso: `PetscFOpen()`, `PetscFClose()`, `PetscPClose()`, `PetscPOpenSetMachine()`

# External Links
$(_doc_external("Sys/PetscPOpen"))
"""
function PetscPOpen(petsclib::PetscLibType, comm::MPI_Comm, machine::Vector{Cchar}, program::Vector{Cchar}, mode::Vector{Cchar}, fp::Libc.FILE) end

@for_petsc function PetscPOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, machine::Vector{Cchar}, program::Vector{Cchar}, mode::Vector{Cchar}, fp::Libc.FILE )

    @chk ccall(
               (:PetscPOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Libc.FILE),
               comm, machine, program, mode, fp,
              )


	return nothing
end 

"""
	PetscPOpenSetMachine(petsclib::PetscLibType,machine::Vector{Cchar}) 
Sets the name of the default machine to run `PetscPOpen()` calls on

Logically Collective, but only the MPI process with rank 0 runs the command

Input Parameter:
- `machine` - machine to run command on or `NULL` for the current machine

Options Database Key:
- `-popen_machine <machine>` - run the process on this machine

Level: intermediate

-seealso: `PetscFOpen()`, `PetscFClose()`, `PetscPClose()`, `PetscPOpen()`

# External Links
$(_doc_external("Sys/PetscPOpenSetMachine"))
"""
function PetscPOpenSetMachine(petsclib::PetscLibType, machine::Vector{Cchar}) end

@for_petsc function PetscPOpenSetMachine(petsclib::$UnionPetscLib, machine::Vector{Cchar} )

    @chk ccall(
               (:PetscPOpenSetMachine, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               machine,
              )


	return nothing
end 

"""
	PetscGetWorkingDirectory(petsclib::PetscLibType,path::Vector{Cchar}, len::Csize_t) 
Gets the current working directory.

Not Collective

Input Parameter:
- `len` - maximum length of `path`

Output Parameter:
- `path` - holds the result value. The string should be long enough to hold the path, for example, `PETSC_MAX_PATH_LEN`

Level: developer

-seealso: `PetscGetTmp()`, `PetscSharedTmp()`, `PetscSharedWorkingDirectory()`, `PetscGetHomeDirectory()`

# External Links
$(_doc_external("Sys/PetscGetWorkingDirectory"))
"""
function PetscGetWorkingDirectory(petsclib::PetscLibType, path::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscGetWorkingDirectory(petsclib::$UnionPetscLib, path::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscGetWorkingDirectory, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               path, len,
              )


	return nothing
end 

"""
	PetscGetHomeDirectory(petsclib::PetscLibType,dir::Vector{Cchar}, maxlen::Csize_t) 
Returns the name of the user's home directory

Not Collective

Input Parameter:
- `maxlen` - maximum length allowed

Output Parameter:
- `dir` - contains the home directory. Must be long enough to hold the name.

Level: developer

-seealso: `PetscGetTmp()`, `PetscSharedTmp()`, `PetscGetWorkingDirectory()`

# External Links
$(_doc_external("Sys/PetscGetHomeDirectory"))
"""
function PetscGetHomeDirectory(petsclib::PetscLibType, dir::Vector{Cchar}, maxlen::Csize_t) end

@for_petsc function PetscGetHomeDirectory(petsclib::$UnionPetscLib, dir::Vector{Cchar}, maxlen::Csize_t )

    @chk ccall(
               (:PetscGetHomeDirectory, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               dir, maxlen,
              )


	return nothing
end 

"""
	PetscFixFilename(petsclib::PetscLibType,filein::Vector{Cchar}, fileout::Vector{Cchar}) 
Fixes a file name so that it is correct for both Unix and
Microsoft Windows by using the correct `/` or `\` to separate directories.

Not Collective

Input Parameter:
- `filein` - name of file to be fixed

Output Parameter:
- `fileout` - the fixed name. Should long enough to hold the filename.

Level: developer

-seealso: `PetscFOpen()`

# External Links
$(_doc_external("Sys/PetscFixFilename"))
"""
function PetscFixFilename(petsclib::PetscLibType, filein::Vector{Cchar}, fileout::Vector{Cchar}) end

@for_petsc function PetscFixFilename(petsclib::$UnionPetscLib, filein::Vector{Cchar}, fileout::Vector{Cchar} )

    @chk ccall(
               (:PetscFixFilename, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               filein, fileout,
              )


	return nothing
end 

"""
	PetscGetRelativePath(petsclib::PetscLibType,fullpath::Vector{Cchar}, path::Vector{Cchar}, flen::Csize_t) 
Given a filename, returns the relative path (removes
all directory specifiers).

Not Collective; No Fortran Support

Input Parameters:
- `fullpath` - full pathname
- `flen`     - size of `path`

Output Parameter:
- `path` - buffer that holds relative pathname

Level: developer

-seealso: `PetscGetFullPath()`

# External Links
$(_doc_external("Sys/PetscGetRelativePath"))
"""
function PetscGetRelativePath(petsclib::PetscLibType, fullpath::Vector{Cchar}, path::Vector{Cchar}, flen::Csize_t) end

@for_petsc function PetscGetRelativePath(petsclib::$UnionPetscLib, fullpath::Vector{Cchar}, path::Vector{Cchar}, flen::Csize_t )

    @chk ccall(
               (:PetscGetRelativePath, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               fullpath, path, flen,
              )


	return nothing
end 

"""
	PetscFormatConvertGetSize(petsclib::PetscLibType,format::Vector{Cchar}, size::Csize_t) 
Gets the length of a string needed to hold data converted with `PetscFormatConvert()` based on the format

No Fortran Support

Input Parameter:
- `format` - the PETSc format string

Output Parameter:
- `size` - the needed length of the new format

Level: developer

-seealso: `PetscFormatConvert()`, `PetscVSNPrintf()`, `PetscVFPrintf()`

# External Links
$(_doc_external("Sys/PetscFormatConvertGetSize"))
"""
function PetscFormatConvertGetSize(petsclib::PetscLibType, format::Vector{Cchar}, size::Csize_t) end

@for_petsc function PetscFormatConvertGetSize(petsclib::$UnionPetscLib, format::Vector{Cchar}, size::Csize_t )

    @chk ccall(
               (:PetscFormatConvertGetSize, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Csize_t}),
               format, size,
              )


	return nothing
end 

"""
	PetscFormatConvert(petsclib::PetscLibType,format::Vector{Cchar}, newformat::Vector{Cchar}) 
converts %g to [|%g|] so that `PetscVSNPrintf()` can ensure all %g formatted numbers have a decimal point when printed.

No Fortran Support

Input Parameter:
- `format` - the PETSc format string

Output Parameter:
- `newformat` - the formatted string, must be long enough to hold result

Level: developer

-seealso: `PetscFormatConvertGetSize()`, `PetscVSNPrintf()`, `PetscVFPrintf()`

# External Links
$(_doc_external("Sys/PetscFormatConvert"))
"""
function PetscFormatConvert(petsclib::PetscLibType, format::Vector{Cchar}, newformat::Vector{Cchar}) end

@for_petsc function PetscFormatConvert(petsclib::$UnionPetscLib, format::Vector{Cchar}, newformat::Vector{Cchar} )

    @chk ccall(
               (:PetscFormatConvert, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               format, newformat,
              )


	return nothing
end 

"""
	PetscFFlush(petsclib::PetscLibType,fd::Libc.FILE) 
Flush a file stream

Input Parameter:
- `fd` - The file stream handle

Level: intermediate

-seealso: `PetscPrintf()`, `PetscFPrintf()`, `PetscVFPrintf()`, `PetscVSNPrintf()`

# External Links
$(_doc_external("Sys/PetscFFlush"))
"""
function PetscFFlush(petsclib::PetscLibType, fd::Libc.FILE) end

@for_petsc function PetscFFlush(petsclib::$UnionPetscLib, fd::Libc.FILE )

    @chk ccall(
               (:PetscFFlush, $petsc_library),
               PetscErrorCode,
               (Ptr{Libc.FILE},),
               fd,
              )


	return nothing
end 

"""
	PetscSynchronizedFlush(petsclib::PetscLibType,comm::MPI_Comm, fd::Libc.FILE) 
Flushes to the screen output from all processors
involved in previous `PetscSynchronizedPrintf()`/`PetscSynchronizedFPrintf()` calls.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `fd`   - the file pointer (valid on MPI rank 0 of the communicator), `PETSC_STDOUT` or value obtained from `PetscFOpen()`

Level: intermediate

-seealso: `PetscSynchronizedPrintf()`, `PetscFPrintf()`, `PetscPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIISynchronizedPrintf()`

# External Links
$(_doc_external("Sys/PetscSynchronizedFlush"))
"""
function PetscSynchronizedFlush(petsclib::PetscLibType, comm::MPI_Comm, fd::Libc.FILE) end

@for_petsc function PetscSynchronizedFlush(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Libc.FILE )

    @chk ccall(
               (:PetscSynchronizedFlush, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}),
               comm, fd,
              )


	return nothing
end 

"""
	PetscSynchronizedFGets(petsclib::PetscLibType,comm::MPI_Comm, fp::Libc.FILE, len::Csize_t, string::Vector{Cchar}) 
Multiple MPI processes all get the same line from a file.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `fp`   - the file pointer
- `len`  - the length of `string`

Output Parameter:
- `string` - the line read from the file, at end of file `string`[0] == 0

Level: intermediate

-seealso: `PetscSynchronizedPrintf()`, `PetscSynchronizedFlush()`,
`PetscFOpen()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIPrintf()`

# External Links
$(_doc_external("Sys/PetscSynchronizedFGets"))
"""
function PetscSynchronizedFGets(petsclib::PetscLibType, comm::MPI_Comm, fp::Libc.FILE, len::Csize_t, string::Vector{Cchar}) end

@for_petsc function PetscSynchronizedFGets(petsclib::$UnionPetscLib, comm::MPI_Comm, fp::Libc.FILE, len::Csize_t, string::Vector{Cchar} )

    @chk ccall(
               (:PetscSynchronizedFGets, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}, Csize_t, Ptr{Cchar}),
               comm, fp, len, string,
              )


	return nothing
end 

"""
	PetscFormatRealArray(petsclib::PetscLibType,buf::Vector{Cchar}, len::Csize_t, fmt::Vector{Cchar}, n::PetscInt, x::Vector{PetscReal}) 

# External Links
$(_doc_external("Sys/PetscFormatRealArray"))
"""
function PetscFormatRealArray(petsclib::PetscLibType, buf::Vector{Cchar}, len::Csize_t, fmt::Vector{Cchar}, n::PetscInt, x::Vector{PetscReal}) end

@for_petsc function PetscFormatRealArray(petsclib::$UnionPetscLib, buf::Vector{Cchar}, len::Csize_t, fmt::Vector{Cchar}, n::$PetscInt, x::Vector{$PetscReal} )

    @chk ccall(
               (:PetscFormatRealArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t, Ptr{Cchar}, $PetscInt, Ptr{$PetscReal}),
               buf, len, fmt, n, x,
              )


	return nothing
end 

"""
	PetscMkdir(petsclib::PetscLibType,dir::Vector{Cchar}) 
Create a directory

Not Collective

Input Parameter:
- `dir` - the directory name

Level: advanced

-seealso: `PetscMkdtemp()`, `PetscRMTree()`

# External Links
$(_doc_external("Sys/PetscMkdir"))
"""
function PetscMkdir(petsclib::PetscLibType, dir::Vector{Cchar}) end

@for_petsc function PetscMkdir(petsclib::$UnionPetscLib, dir::Vector{Cchar} )

    @chk ccall(
               (:PetscMkdir, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               dir,
              )


	return nothing
end 

"""
	PetscMkdtemp(petsclib::PetscLibType,dir::Vector{Cchar}) 
Create a directory with a unique name given a name template.

Input Parameter:
- `dir` - file name template, the last six characters must be 'XXXXXX', and they will be modified upon return

Level: advanced

-seealso: `PetscMkdir()`, `PetscRMTree()`

# External Links
$(_doc_external("Sys/PetscMkdtemp"))
"""
function PetscMkdtemp(petsclib::PetscLibType, dir::Vector{Cchar}) end

@for_petsc function PetscMkdtemp(petsclib::$UnionPetscLib, dir::Vector{Cchar} )

    @chk ccall(
               (:PetscMkdtemp, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               dir,
              )


	return nothing
end 

"""
	PetscRMTree(petsclib::PetscLibType,dir::Vector{Cchar}) 
delete a directory and all of its children

Input Parameter:
- `dir` - the name of the directory

Level: advanced

-seealso: `PetscMkdtemp()`, `PetscMkdir()`

# External Links
$(_doc_external("Sys/PetscRMTree"))
"""
function PetscRMTree(petsclib::PetscLibType, dir::Vector{Cchar}) end

@for_petsc function PetscRMTree(petsclib::$UnionPetscLib, dir::Vector{Cchar} )

    @chk ccall(
               (:PetscRMTree, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               dir,
              )


	return nothing
end 

"""
	PetscPythonFinalize(petsclib::PetscLibType) 
Finalize PETSc for use with Python.

Level: intermediate

-seealso: `PetscPythonInitialize()`, `PetscPythonPrintError()`

# External Links
$(_doc_external("Sys/PetscPythonFinalize"))
"""
function PetscPythonFinalize(petsclib::PetscLibType) end

@for_petsc function PetscPythonFinalize(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPythonFinalize, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscPythonInitialize(petsclib::PetscLibType,pyexe::Vector{Cchar}, pylib::Vector{Cchar}) 
Initialize Python for use with PETSc and import petsc4py.

Input Parameters:
- `pyexe`  - path to the Python interpreter executable, or `NULL`.
- `pylib`  - full path to the Python dynamic library, or `NULL`.

Options Database Key:
- `-python <exe>` - Initializes Python, and optionally takes a Python executable name

Level: intermediate

-seealso: `PetscPythonFinalize()`, `PetscPythonPrintError()`

# External Links
$(_doc_external("Sys/PetscPythonInitialize"))
"""
function PetscPythonInitialize(petsclib::PetscLibType, pyexe::Vector{Cchar}, pylib::Vector{Cchar}) end

@for_petsc function PetscPythonInitialize(petsclib::$UnionPetscLib, pyexe::Vector{Cchar}, pylib::Vector{Cchar} )

    @chk ccall(
               (:PetscPythonInitialize, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               pyexe, pylib,
              )


	return nothing
end 

"""
	PetscPythonPrintError(petsclib::PetscLibType) 
Print any current Python errors.

Level: developer

-seealso: `PetscPythonInitialize()`, `PetscPythonFinalize()`

# External Links
$(_doc_external("Sys/PetscPythonPrintError"))
"""
function PetscPythonPrintError(petsclib::PetscLibType) end

@for_petsc function PetscPythonPrintError(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPythonPrintError, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscPythonMonitorSet(petsclib::PetscLibType,obj::PetscObject, url::Vector{Cchar}) 
Set a Python monitor for a `PetscObject`

Level: developer

-seealso: `PetscPythonInitialize()`, `PetscPythonFinalize()`, `PetscPythonPrintError()`

# External Links
$(_doc_external("Sys/PetscPythonMonitorSet"))
"""
function PetscPythonMonitorSet(petsclib::PetscLibType, obj::PetscObject, url::Vector{Cchar}) end

@for_petsc function PetscPythonMonitorSet(petsclib::$UnionPetscLib, obj::PetscObject, url::Vector{Cchar} )

    @chk ccall(
               (:PetscPythonMonitorSet, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, url,
              )


	return nothing
end 

"""
	PetscSysFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the system library portion of PETSc.
It is called from `PetscFinalize()`.

Level: developer

-seealso: `PetscSysInitializePackage()`, `PetscFinalize()`

# External Links
$(_doc_external("Sys/PetscSysFinalizePackage"))
"""
function PetscSysFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscSysFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSysFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscSysInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the system library portion of PETSc. It is called
from `PetscDLLibraryRegister_petsc()` when using dynamic libraries, and in the call to `PetscInitialize()`
when using shared or static libraries.

Level: developer

-seealso: `PetscSysFinalizePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscSysInitializePackage"))
"""
function PetscSysInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscSysInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSysInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	identical::PetscBool = PetscMonitorCompare(petsclib::PetscLibType,nmon::external, nmctx::Cvoid, nmdestroy::PetscCtxDestroyFn, mon::external, mctx::Cvoid, mdestroy::PetscCtxDestroyFn) 
Checks if two monitors are identical; if they are then it destroys the new one

Not Collective

Input Parameters:
- `nmon`      - The new monitor
- `nmctx`     - The new monitor context, or `NULL`
- `nmdestroy` - The new monitor context destroy function, or `NULL`, see `PetscCtxDestroyFn` for its calling sequence
- `mon`       - The old monitor
- `mctx`      - The old monitor context, or `NULL`
- `mdestroy`  - The old monitor context destroy function, or `NULL`, see `PetscCtxDestroyFn` for its calling sequence

Output Parameter:
- `identical` - `PETSC_TRUE` if the monitors are the same

Level: developer

-seealso: [](sec_viewers), `DMMonitorSetFromOptions()`, `KSPMonitorSetFromOptions()`, `SNESMonitorSetFromOptions()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Sys/PetscMonitorCompare"))
"""
function PetscMonitorCompare(petsclib::PetscLibType, nmon::external, nmctx::Cvoid, nmdestroy::PetscCtxDestroyFn, mon::external, mctx::Cvoid, mdestroy::PetscCtxDestroyFn) end

@for_petsc function PetscMonitorCompare(petsclib::$UnionPetscLib, nmon::external, nmctx::Cvoid, nmdestroy::PetscCtxDestroyFn, mon::external, mctx::Cvoid, mdestroy::PetscCtxDestroyFn )
	identical_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscMonitorCompare, $petsc_library),
               PetscErrorCode,
               (external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}, Ptr{PetscBool}),
               nmon, nmctx, nmdestroy, mon, mctx, mdestroy, identical_,
              )

	identical = identical_[]

	return identical
end 

"""
	htype::hid_t = PetscDataTypeToHDF5DataType(petsclib::PetscLibType,ptype::PetscDataType) 
Converts the PETSc name of a datatype to its HDF5 name.

Not Collective

Input Parameter:
- `ptype` - the PETSc datatype name (for example `PETSC_DOUBLE`)

Output Parameter:
- `htype` - the HDF5  datatype

Level: advanced

-seealso: [](sec_viewers), `PetscDataType`, `PetscHDF5DataTypeToPetscDataType()`

# External Links
$(_doc_external("Sys/PetscDataTypeToHDF5DataType"))
"""
function PetscDataTypeToHDF5DataType(petsclib::PetscLibType, ptype::PetscDataType) end

@for_petsc function PetscDataTypeToHDF5DataType(petsclib::$UnionPetscLib, ptype::PetscDataType )
	htype_ = Ref{hid_t}()

    @chk ccall(
               (:PetscDataTypeToHDF5DataType, $petsc_library),
               PetscErrorCode,
               (PetscDataType, Ptr{hid_t}),
               ptype, htype_,
              )

	htype = htype_[]

	return htype
end 

"""
	ptype::PetscDataType = PetscHDF5DataTypeToPetscDataType(petsclib::PetscLibType,htype::hid_t) 
Finds the PETSc name of a datatype from its HDF5 name

Not Collective

Input Parameter:
- `htype` - the HDF5 datatype (for example `H5T_NATIVE_DOUBLE`, ...)

Output Parameter:
- `ptype` - the PETSc datatype name (for example `PETSC_DOUBLE`)

Level: advanced

-seealso: [](sec_viewers), `PetscDataType`

# External Links
$(_doc_external("Sys/PetscHDF5DataTypeToPetscDataType"))
"""
function PetscHDF5DataTypeToPetscDataType(petsclib::PetscLibType, htype::hid_t) end

@for_petsc function PetscHDF5DataTypeToPetscDataType(petsclib::$UnionPetscLib, htype::hid_t )
	ptype_ = Ref{PetscDataType}()

    @chk ccall(
               (:PetscHDF5DataTypeToPetscDataType, $petsc_library),
               PetscErrorCode,
               (hid_t, Ptr{PetscDataType}),
               htype, ptype_,
              )

	ptype = unsafe_string(ptype_[])

	return ptype
end 

"""
	PetscOpenSocket(petsclib::PetscLibType,hostname::Vector{Cchar}, portnum::Cint, t::Cint) 
handles connected to an open port where someone is waiting.

Input Parameters:
- `hostname` - for example www.mcs.anl.gov
- `portnum`  - for example 80

Output Parameter:
- `t` - the socket number

-seealso: `PetscSocketListen()`, `PetscSocketEstablish()`, `PetscHTTPRequest()`, `PetscHTTPSConnect()`

# External Links
$(_doc_external("Sys/PetscOpenSocket"))
"""
function PetscOpenSocket(petsclib::PetscLibType, hostname::Vector{Cchar}, portnum::Cint, t::Cint) end

@for_petsc function PetscOpenSocket(petsclib::$UnionPetscLib, hostname::Vector{Cchar}, portnum::Cint, t::Cint )

    @chk ccall(
               (:PetscOpenSocket, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cint, Ptr{Cint}),
               hostname, portnum, t,
              )


	return nothing
end 

"""
	PetscBTView(petsclib::PetscLibType,m::PetscCount, bt::PetscBT, viewer::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscBTView"))
"""
function PetscBTView(petsclib::PetscLibType, m::PetscCount, bt::PetscBT, viewer::PetscViewer) end

@for_petsc function PetscBTView(petsclib::$UnionPetscLib, m::PetscCount, bt::PetscBT, viewer::PetscViewer )

    @chk ccall(
               (:PetscBTView, $petsc_library),
               PetscErrorCode,
               (PetscCount, PetscBT, PetscViewer),
               m, bt, viewer,
              )


	return nothing
end 

"""
	PetscRegisterFinalize(petsclib::PetscLibType,f::external) 
Registers a function that is to be called in `PetscFinalize()`

Not Collective

Input Parameter:
- `f` - function to be called

Level: developer

-seealso: `PetscRegisterFinalizeAll()`, `PetscObjectRegisterDestroy()`

# External Links
$(_doc_external("Sys/PetscRegisterFinalize"))
"""
function PetscRegisterFinalize(petsclib::PetscLibType, f::external) end

@for_petsc function PetscRegisterFinalize(petsclib::$UnionPetscLib, f::external )

    @chk ccall(
               (:PetscRegisterFinalize, $petsc_library),
               PetscErrorCode,
               (external,),
               f,
              )


	return nothing
end 

"""
	PetscRegisterFinalizeAll(petsclib::PetscLibType) 
Runs all the finalize functions set with `PetscRegisterFinalize()`

Not Collective unless registered functions are collective

Level: developer

-seealso: `PetscRegisterFinalize()`, `PetscObjectRegisterDestroyAll()`

# External Links
$(_doc_external("Sys/PetscRegisterFinalizeAll"))
"""
function PetscRegisterFinalizeAll(petsclib::PetscLibType) end

@for_petsc function PetscRegisterFinalizeAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscRegisterFinalizeAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	has::PetscBool = PetscHasExternalPackage(petsclib::PetscLibType,pkg::Vector{Cchar}) 
Determine whether PETSc has been configured with the given package

Not Collective

Input Parameter:
- `pkg` - external package name

Output Parameter:
- `has` - `PETSC_TRUE` if PETSc is configured with the given package, else `PETSC_FALSE`.

Level: intermediate

-seealso: `PetscViewerType`, `MatPartitioningType`, `MatSolverType`

# External Links
$(_doc_external("Sys/PetscHasExternalPackage"))
"""
function PetscHasExternalPackage(petsclib::PetscLibType, pkg::Vector{Cchar}) end

@for_petsc function PetscHasExternalPackage(petsclib::$UnionPetscLib, pkg::Vector{Cchar} )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscHasExternalPackage, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscBool}),
               pkg, has_,
              )

	has = has_[]

	return has
end 

"""
	PetscInitializeNoPointers(petsclib::PetscLibType,argc::Cint, args::Cchar, filename::Vector{Cchar}, help::Vector{Cchar}) 
Calls PetscInitialize() from C/C++ without the pointers to argc and args

Collective, No Fortran Support

Input Parameters:
- `argc`     - number of args
- `args`     - array of command line arguments
- `filename` - optional name of the program file, pass `NULL` to ignore
- `help`     - optional help, pass `NULL` to ignore

Level: advanced

-seealso: `PetscInitialize()`, `PetscInitializeFortran()`, `PetscInitializeNoArguments()`
*/
PetscErrorCode PetscInitializeNoPointers(int argc, char **args, const char *filename, const char *help)
{
int    myargc = argc;
char **myargs = args;

PetscFunctionBegin;
PetscCall(PetscInitialize(&myargc, &myargs, filename, help));
PetscCall(PetscPopSignalHandler());
PetscBeganMPI = PETSC_FALSE;
PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
PetscInitializeNoArguments - Calls `PetscInitialize()` from C/C++ without
the command line arguments.

Collective

Level: advanced

-seealso: `PetscInitialize()`, `PetscInitializeFortran()`

# External Links
$(_doc_external("Sys/PetscInitializeNoPointers"))
"""
function PetscInitializeNoPointers(petsclib::PetscLibType, argc::Cint, args::Cchar, filename::Vector{Cchar}, help::Vector{Cchar}) end

@for_petsc function PetscInitializeNoPointers(petsclib::$UnionPetscLib, argc::Cint, args::Cchar, filename::Vector{Cchar}, help::Vector{Cchar} )

    @chk ccall(
               (:PetscInitializeNoPointers, $petsc_library),
               PetscErrorCode,
               (Cint, Cchar, Ptr{Cchar}, Ptr{Cchar}),
               argc, args, filename, help,
              )


	return nothing
end 

"""
	PetscInitializeNoArguments(petsclib::PetscLibType) 
Calls `PetscInitialize()` from C/C++ without
the command line arguments.

Collective

Level: advanced

-seealso: `PetscInitialize()`, `PetscInitializeFortran()`

# External Links
$(_doc_external("Sys/PetscInitializeNoArguments"))
"""
function PetscInitializeNoArguments(petsclib::PetscLibType) end

@for_petsc function PetscInitializeNoArguments(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscInitializeNoArguments, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	isInitialized::PetscBool = PetscInitialized(petsclib::PetscLibType) 
Determine whether PETSc is initialized.

Output Parameter:
- `isInitialized` - `PETSC_TRUE` if PETSc is initialized, `PETSC_FALSE` otherwise

Level: beginner

-seealso: `PetscInitialize()`, `PetscInitializeNoArguments()`, `PetscInitializeFortran()`

# External Links
$(_doc_external("Sys/PetscInitialized"))
"""
function PetscInitialized(petsclib::PetscLibType) end

@for_petsc function PetscInitialized(petsclib::$UnionPetscLib)
	isInitialized_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscInitialized, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               isInitialized_,
              )

	isInitialized = isInitialized_[]

	return isInitialized
end 

"""
	isFinalized::PetscBool = PetscFinalized(petsclib::PetscLibType) 
Determine whether `PetscFinalize()` has been called yet

Output Parameter:
- `isFinalized` - `PETSC_TRUE` if PETSc is finalized, `PETSC_FALSE` otherwise

Level: developer

-seealso: `PetscInitialize()`, `PetscInitializeNoArguments()`, `PetscInitializeFortran()`

# External Links
$(_doc_external("Sys/PetscFinalized"))
"""
function PetscFinalized(petsclib::PetscLibType) end

@for_petsc function PetscFinalized(petsclib::$UnionPetscLib)
	isFinalized_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscFinalized, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               isFinalized_,
              )

	isFinalized = isFinalized_[]

	return isFinalized
end 

"""
	max::PetscInt,sum::PetscInt = PetscMaxSum(petsclib::PetscLibType,comm::MPI_Comm, array::Vector{PetscInt}) 
Returns the max of the first entry over all MPI processes and the sum of the second entry.

Collective

Input Parameters:
- `comm`  - the communicator
- `array` - an arry of length 2 times `size`, the number of MPI processes

Output Parameters:
- `max` - the maximum of `array[2*rank]` over all MPI processes
- `sum` - the sum of the `array[2*rank + 1]` over all MPI processes

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscMaxSum"))
"""
function PetscMaxSum(petsclib::PetscLibType, comm::MPI_Comm, array::Vector{PetscInt}) end

@for_petsc function PetscMaxSum(petsclib::$UnionPetscLib, comm::MPI_Comm, array::Vector{$PetscInt} )
	max_ = Ref{$PetscInt}()
	sum_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscMaxSum, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               comm, array, max_, sum_,
              )

	max = max_[]
	sum = sum_[]

	return max,sum
end 

"""
	PetscSetProgramName(petsclib::PetscLibType,name::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscSetProgramName"))
"""
function PetscSetProgramName(petsclib::PetscLibType, name::Vector{Cchar}) end

@for_petsc function PetscSetProgramName(petsclib::$UnionPetscLib, name::Vector{Cchar} )

    @chk ccall(
               (:PetscSetProgramName, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               name,
              )


	return nothing
end 

"""
	PetscGetProgramName(petsclib::PetscLibType,name::Vector{Cchar}, len::Csize_t) 
Gets the name of the running program.

Not Collective

Input Parameter:
- `len` - length of the string name

Output Parameter:
- `name` - the name of the running program, provide a string of length `PETSC_MAX_PATH_LEN`

Level: advanced

-seealso: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArguments()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscGetProgramName"))
"""
function PetscGetProgramName(petsclib::PetscLibType, name::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscGetProgramName(petsclib::$UnionPetscLib, name::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscGetProgramName, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               name, len,
              )


	return nothing
end 

"""
	PetscGetArgs(petsclib::PetscLibType,argc::Cint, args::Cchar) 
Allows you to access the raw command line arguments anywhere
after `PetscInitialize()` is called but before `PetscFinalize()`.

Not Collective, No Fortran Support

Output Parameters:
- `argc` - count of the number of command line arguments
- `args` - the command line arguments

Level: intermediate

-seealso: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArguments()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscGetArgs"))
"""
function PetscGetArgs(petsclib::PetscLibType, argc::Cint, args::Cchar) end

@for_petsc function PetscGetArgs(petsclib::$UnionPetscLib, argc::Cint, args::Cchar )

    @chk ccall(
               (:PetscGetArgs, $petsc_library),
               PetscErrorCode,
               (Ptr{Cint}, Cchar),
               argc, args,
              )


	return nothing
end 

"""
	PetscGetArguments(petsclib::PetscLibType,args::Cchar) 
Allows you to access the command line arguments anywhere
after `PetscInitialize()` is called but before `PetscFinalize()`.

Not Collective, No Fortran Support

Output Parameter:
- `args` - the command line arguments

Level: intermediate

-seealso: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArgs()`, `PetscFreeArguments()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscGetArguments"))
"""
function PetscGetArguments(petsclib::PetscLibType, args::Cchar) end

@for_petsc function PetscGetArguments(petsclib::$UnionPetscLib, args::Cchar )

    @chk ccall(
               (:PetscGetArguments, $petsc_library),
               PetscErrorCode,
               (Cchar,),
               args,
              )


	return nothing
end 

"""
	PetscFreeArguments(petsclib::PetscLibType,args::Cchar) 
Frees the memory obtained with `PetscGetArguments()`

Not Collective, No Fortran Support

Output Parameter:
- `args` - the command line arguments

Level: intermediate

-seealso: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArgs()`, `PetscGetArguments()`

# External Links
$(_doc_external("Sys/PetscFreeArguments"))
"""
function PetscFreeArguments(petsclib::PetscLibType, args::Cchar) end

@for_petsc function PetscFreeArguments(petsclib::$UnionPetscLib, args::Cchar )

    @chk ccall(
               (:PetscFreeArguments, $petsc_library),
               PetscErrorCode,
               (Cchar,),
               args,
              )


	return nothing
end 

"""
	PetscInitialize(petsclib::PetscLibType,argc::Cint, args::Cchar, file::Vector{Cchar}, help::Vector{Cchar}) 
Initializes the PETSc database and MPI.
`PetscInitialize()` calls MPI_Init() if that has yet to be called,
so this routine should always be called near the beginning of
your program -- usually the very first line!

Collective on `MPI_COMM_WORLD` or `PETSC_COMM_WORLD` if it has been set

Input Parameters:
- `argc` - count of number of command line arguments
- `args` - the command line arguments
- `file` - [optional] PETSc database file, append ":yaml" to filename to specify YAML options format.
Use `NULL` or empty string to not check for code specific file.
Also checks `~/.petscrc`, `.petscrc` and `petscrc`.
Use `-skip_petscrc` in the code specific file (or command line) to skip `~/.petscrc`, `.petscrc` and `petscrc` files.
- `help` - [optional] Help message to print, use `NULL` for no message

If you wish PETSc code to run ONLY on a subcommunicator of `MPI_COMM_WORLD`, create that
communicator first and assign it to `PETSC_COMM_WORLD` BEFORE calling `PetscInitialize()`.
then do this. If ALL processes in the job are using `PetscInitialize()` and `PetscFinalize()` then you don't need to do this, even
if different subcommunicators of the job are doing different things with PETSc.

Options Database Keys:
- `-help [intro]`                                       - prints help method for each option; if `intro` is given the program stops after printing the introductory help message
- `-start_in_debugger [noxterm,dbx,xdb,gdb,...]`        - Starts program in debugger
- `-on_error_attach_debugger [noxterm,dbx,xdb,gdb,...]` - Starts debugger when error detected
- `-on_error_emacs <machinename>`                       - causes `emacsclient` to jump to error file if an error is detected
- `-on_error_abort`                                     - calls `abort()` when error detected (no traceback)
- `-on_error_mpiabort`                                  - calls `MPI_abort()` when error detected
- `-error_output_stdout`                                - prints PETSc error messages to `stdout` instead of the default `stderr`
- `-error_output_none`                                  - does not print the error messages (but handles errors in the same way as if this was not called)
- `-debugger_ranks [rank1,rank2,...]`                   - Indicates MPI ranks to start in debugger
- `-debugger_pause [sleeptime] (in seconds)`            - Pauses debugger, use if it takes a long time for the debugger to start up on your system
- `-stop_for_debugger`                                  - Print message on how to attach debugger manually to
process and wait (`-debugger_pause`) seconds for attachment
- `-malloc_dump`                                        - prints a list of all unfreed memory at the end of the run
- `-malloc_test`                                        - like `-malloc_dump` `-malloc_debug`, only active for debugging build, ignored in optimized build. Often set in `PETSC_OPTIONS` environmental variable
- `-malloc_view`                                        - show a list of all allocated memory during `PetscFinalize()`
- `-malloc_view_threshold <t>`                          - only list memory allocations of size greater than t with `-malloc_view`
- `-malloc_requested_size`                              - malloc logging will record the requested size rather than (possibly large) size after alignment
- `-fp_trap`                                            - Stops on floating point exceptions
- `-no_signal_handler`                                  - Indicates not to trap error signals
- `-shared_tmp`                                         - indicates `/tmp` directory is known to be shared by all processors
- `-not_shared_tmp`                                     - indicates each processor has own `/tmp`
- `-tmp`                                                - alternative directory to use instead of `/tmp`
- `-python <exe>`                                       - Initializes Python, and optionally takes a Python executable name
- `-mpiuni-allow-multiprocess-launch`                   - allow `mpiexec` to launch multiple independent MPI-Uni jobs, otherwise a sanity check error is invoked to prevent misuse of MPI-Uni

Options Database Keys for Option Database:
- `-skip_petscrc`           - skip the default option files `~/.petscrc`, `.petscrc`, `petscrc`
- `-options_monitor`        - monitor all set options to standard output for the whole program run
- `-options_monitor_cancel` - cancel options monitoring hard-wired using `PetscOptionsMonitorSet()`

Options -options_monitor_{all,cancel} are
position-independent and apply to all options set since the PETSc start.
They can be used also in option files.

See `PetscOptionsMonitorSet()` to do monitoring programmatically.

Options Database Keys for Profiling:
See Users-Manual: ch_profiling for details.
- `-info [filename][:[~]<list,of,classnames>[:[~]self]]` - Prints verbose information. See `PetscInfo()`.
- `-log_sync`                                            - Enable barrier synchronization for all events. This option is useful to debug imbalance within each event,
however it slows things down and gives a distorted view of the overall runtime.
- `-log_trace [filename]`                                - Print traces of all PETSc calls to the screen (useful to determine where a program
hangs without running in the debugger).  See `PetscLogTraceBegin()`.
- `-log_view [:filename:format][,[:filename:format]...]` - Prints summary of flop and timing information to screen or file, see `PetscLogView()` (up to 4 viewers)
- `-log_view_memory`                                     - Includes in the summary from -log_view the memory used in each event, see `PetscLogView()`.
- `-log_view_gpu_time`                                   - Includes in the summary from -log_view the time used in each GPU kernel, see `PetscLogView().
- `-log_exclude: <vec,mat,pc,ksp,snes>`                  - excludes subset of object classes from logging
- `-log [filename]`                                      - Logs profiling information in a dump file, see `PetscLogDump()`.
- `-log_all [filename]`                                  - Same as `-log`.
- `-log_mpe [filename]`                                  - Creates a logfile viewable by the utility Jumpshot (in MPICH distribution)
- `-log_perfstubs`                                       - Starts a log handler with the perfstubs interface (which is used by TAU)
- `-log_nvtx`                                            - Starts an nvtx log handler for use with Nsight
- `-log_roctx`                                           - Starts an roctx log handler for use with rocprof on AMD GPUs
- `-viewfromoptions on,off`                              - Enable or disable `XXXSetFromOptions()` calls, for applications with many small solves turn this off
- `-get_total_flops`                                     - Returns total flops done by all processors
- `-memory_view`                                         - Print memory usage at end of run
- `-check_pointer_intensity 0,1,2`                       - if pointers are checked for validity (debug version only), using 0 will result in faster code

Options Database Keys for SAWs:
- `-saws_port <portnumber>`        - port number to publish SAWs data, default is 8080
- `-saws_port_auto_select`         - have SAWs select a new unique port number where it publishes the data, the URL is printed to the screen
this is useful when you are running many jobs that utilize SAWs at the same time
- `-saws_log <filename>`           - save a log of all SAWs communication
- `-saws_https <certificate file>` - have SAWs use HTTPS instead of HTTP
- `-saws_root <directory>`         - allow SAWs to have access to the given directory to search for requested resources and files

Environmental Variables:
- `PETSC_TMP`                     - alternative directory to use instead of `/tmp`
- `PETSC_SHARED_TMP`              - `/tmp` is shared by all processes
- `PETSC_NOT_SHARED_TMP`          - each process has its own private `/tmp`
- `PETSC_OPTIONS`                 - a string containing additional options for PETSc in the form of command line "-key value" pairs
- `PETSC_OPTIONS_YAML`            - (requires configuring PETSc to use libyaml with `--download-yaml`) a string containing additional options for PETSc in the form of a YAML document
- `PETSC_VIEWER_SOCKET_PORT`      - socket number to use for socket viewer
- `PETSC_VIEWER_SOCKET_MACHINE`   - machine to use for socket viewer to connect to

Level: beginner

-seealso: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArgs()`, `PetscInitializeNoArguments()`, `PetscLogGpuTime()`

# External Links
$(_doc_external("Sys/PetscInitialize"))
"""
function PetscInitialize(petsclib::PetscLibType, argc::Cint, args::Cchar, file::Vector{Cchar}, help::Vector{Cchar}) end

@for_petsc function PetscInitialize(petsclib::$UnionPetscLib, argc::Cint, args::Cchar, file::Vector{Cchar}, help::Vector{Cchar} )

    @chk ccall(
               (:PetscInitialize, $petsc_library),
               PetscErrorCode,
               (Ptr{Cint}, Cchar, Ptr{Cchar}, Ptr{Cchar}),
               argc, args, file, help,
              )


	return nothing
end 

"""
	PetscFinalize(petsclib::PetscLibType) 
Checks for options to be called at the conclusion of a PETSc program and frees any remaining PETSc objects and data structures.
of the program. Automatically calls `MPI_Finalize()` if the user had not called `MPI_Init()` before calling `PetscInitialize()`.

Collective on `PETSC_COMM_WORLD`

Options Database Keys:
- `-options_view`                    - Calls `PetscOptionsView()`
- `-options_left`                    - Prints unused options that remain in the database
- `-objects_dump [all]`              - Prints list of objects allocated by the user that have not been freed, the option all cause all outstanding objects to be listed
- `-mpidump`                         - Calls PetscMPIDump()
- `-malloc_dump <optional filename>` - Calls `PetscMallocDump()`, displays all memory allocated that has not been freed
- `-memory_view`                     - Prints total memory usage
- `-malloc_view <optional filename>` - Prints list of all memory allocated and in what functions

Level: beginner

-seealso: `PetscInitialize()`, `PetscOptionsView()`, `PetscMallocDump()`, `PetscMPIDump()`, `PetscEnd()`

# External Links
$(_doc_external("Sys/PetscFinalize"))
"""
function PetscFinalize(petsclib::PetscLibType) end

@for_petsc function PetscFinalize(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFinalize, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscEnd(petsclib::PetscLibType) 
Calls `PetscFinalize()` and then ends the program. This is useful if one
wishes a clean exit somewhere deep in the program.

Collective on `PETSC_COMM_WORLD`

Level: advanced

-seealso: `PetscInitialize()`, `PetscOptionsView()`, `PetscMallocDump()`, `PetscMPIDump()`, `PetscFinalize()`

# External Links
$(_doc_external("Sys/PetscEnd"))
"""
function PetscEnd(petsclib::PetscLibType) end

@for_petsc function PetscEnd(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscEnd, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscSetHelpVersionFunctions(petsclib::PetscLibType,help::external, version::external) 
Sets functions that print help and version information
before the PETSc help and version information is printed.

No Fortran Support

Input Parameters:
- `help`    - the help function (may be `NULL`)
- `version` - the version function (may be `NULL`)

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscSetHelpVersionFunctions"))
"""
function PetscSetHelpVersionFunctions(petsclib::PetscLibType, help::external, version::external) end

@for_petsc function PetscSetHelpVersionFunctions(petsclib::$UnionPetscLib, help::external, version::external )

    @chk ccall(
               (:PetscSetHelpVersionFunctions, $petsc_library),
               PetscErrorCode,
               (external, external),
               help, version,
              )


	return nothing
end 

"""
	PetscInitializeFortran(petsclib::PetscLibType) 
Routine that should be called soon AFTER
the call to `PetscInitialize()` if one is using a C main program
that calls Fortran routines that in turn call PETSc routines.

Collective on `PETSC_COMM_WORLD`

Level: beginner

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscInitializeFortran"))
"""
function PetscInitializeFortran(petsclib::PetscLibType) end

@for_petsc function PetscInitializeFortran(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscInitializeFortran, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscCommGetNewTag(petsclib::PetscLibType,comm::MPI_Comm, tag::PetscMPIInt) 
Gets a unique new tag from a PETSc communicator

Collective

Input Parameter:
- `comm` - the MPI communicator

Output Parameter:
- `tag` - the new tag

Level: developer

-seealso: `PetscObjectGetNewTag()`, `PetscCommDuplicate()`

# External Links
$(_doc_external("Sys/PetscCommGetNewTag"))
"""
function PetscCommGetNewTag(petsclib::PetscLibType, comm::MPI_Comm, tag::PetscMPIInt) end

@for_petsc function PetscCommGetNewTag(petsclib::$UnionPetscLib, comm::MPI_Comm, tag::PetscMPIInt )

    @chk ccall(
               (:PetscCommGetNewTag, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscMPIInt}),
               comm, tag,
              )


	return nothing
end 

"""
	PetscCommGetComm(petsclib::PetscLibType,comm_in::MPI_Comm, comm_out::MPI_Comm) 
get a new MPI communicator from a PETSc communicator that can be passed off to another package

Collective

Input Parameter:
- `comm_in` - Input communicator

Output Parameter:
- `comm_out` - Output communicator

Level: developer

-seealso: `PetscObjectGetNewTag()`, `PetscCommGetNewTag()`, `PetscCommDestroy()`, `PetscCommRestoreComm()`

# External Links
$(_doc_external("Sys/PetscCommGetComm"))
"""
function PetscCommGetComm(petsclib::PetscLibType, comm_in::MPI_Comm, comm_out::MPI_Comm) end

@for_petsc function PetscCommGetComm(petsclib::$UnionPetscLib, comm_in::MPI_Comm, comm_out::MPI_Comm )

    @chk ccall(
               (:PetscCommGetComm, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{MPI_Comm}),
               comm_in, comm_out,
              )


	return nothing
end 

"""
	PetscCommRestoreComm(petsclib::PetscLibType,comm_in::MPI_Comm, comm_out::MPI_Comm) 
restores an MPI communicator that was obtained with `PetscCommGetComm()`

Collective

Input Parameters:
- `comm_in`  - Input communicator
- `comm_out` - returned communicator

Level: developer

-seealso: `PetscObjectGetNewTag()`, `PetscCommGetNewTag()`, `PetscCommDestroy()`

# External Links
$(_doc_external("Sys/PetscCommRestoreComm"))
"""
function PetscCommRestoreComm(petsclib::PetscLibType, comm_in::MPI_Comm, comm_out::MPI_Comm) end

@for_petsc function PetscCommRestoreComm(petsclib::$UnionPetscLib, comm_in::MPI_Comm, comm_out::MPI_Comm )

    @chk ccall(
               (:PetscCommRestoreComm, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{MPI_Comm}),
               comm_in, comm_out,
              )


	return nothing
end 

"""
	comm_out::MPI_Comm,first_tag::PetscMPIInt = PetscCommDuplicate(petsclib::PetscLibType,comm_in::MPI_Comm) 
Duplicates the communicator only if it is not already a PETSc communicator.

Collective

Input Parameter:
- `comm_in` - Input communicator

Output Parameters:
- `comm_out`  - Output communicator.  May be `comm_in`.
- `first_tag` - Tag available that has not already been used with this communicator (you may pass in `NULL` if you do not need a tag)

Level: developer

-seealso: `PetscObjectGetNewTag()`, `PetscCommGetNewTag()`, `PetscCommDestroy()`

# External Links
$(_doc_external("Sys/PetscCommDuplicate"))
"""
function PetscCommDuplicate(petsclib::PetscLibType, comm_in::MPI_Comm) end

@for_petsc function PetscCommDuplicate(petsclib::$UnionPetscLib, comm_in::MPI_Comm )
	comm_out_ = Ref{MPI_Comm}()
	first_tag_ = Ref{PetscMPIInt}()

    @chk ccall(
               (:PetscCommDuplicate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{MPI_Comm}, Ptr{PetscMPIInt}),
               comm_in, comm_out_, first_tag_,
              )

	comm_out = comm_out_[]
	first_tag = first_tag_[]

	return comm_out,first_tag
end 

"""
	PetscCommDestroy(petsclib::PetscLibType,comm::MPI_Comm) 
Frees communicator obtained with `PetscCommDuplicate()`.

Collective

Input Parameter:
- `comm` - the communicator to free

Level: developer

-seealso: `PetscCommDuplicate()`

# External Links
$(_doc_external("Sys/PetscCommDestroy"))
"""
function PetscCommDestroy(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscCommDestroy(petsclib::$UnionPetscLib, comm::MPI_Comm )

    @chk ccall(
               (:PetscCommDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{MPI_Comm},),
               comm,
              )


	return nothing
end 

"""
	PetscGetVersion(petsclib::PetscLibType,version::Vector{Cchar}, len::Csize_t) 
Gets the PETSc version information in a string.

Not Collective; No Fortran Support

Input Parameter:
- `len` - length of the string

Output Parameter:
- `version` - version string

Level: developer

-seealso: `PetscGetProgramName()`, `PetscGetVersionNumber()`

# External Links
$(_doc_external("Sys/PetscGetVersion"))
"""
function PetscGetVersion(petsclib::PetscLibType, version::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscGetVersion(petsclib::$UnionPetscLib, version::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscGetVersion, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               version, len,
              )


	return nothing
end 

"""
	major::PetscInt,minor::PetscInt,subminor::PetscInt,release::PetscInt = PetscGetVersionNumber(petsclib::PetscLibType) 
Gets the PETSc version information from the library

Not Collective

Output Parameters:
- `major`    - the major version (optional, pass `NULL` if not requested)
- `minor`    - the minor version (optional, pass `NULL` if not requested)
- `subminor` - the subminor version (patch number)  (optional, pass `NULL` if not requested)
- `release`  - indicates the library is from a release, not random git repository  (optional, pass `NULL` if not requested)

Level: developer

-seealso: `PetscGetProgramName()`, `PetscGetVersion()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscGetVersionNumber"))
"""
function PetscGetVersionNumber(petsclib::PetscLibType) end

@for_petsc function PetscGetVersionNumber(petsclib::$UnionPetscLib)
	major_ = Ref{$PetscInt}()
	minor_ = Ref{$PetscInt}()
	subminor_ = Ref{$PetscInt}()
	release_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscGetVersionNumber, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               major_, minor_, subminor_, release_,
              )

	major = major_[]
	minor = minor_[]
	subminor = subminor_[]
	release = release_[]

	return major,minor,subminor,release
end 

"""
	PetscBLASSetNumThreads(petsclib::PetscLibType,nt::PetscInt) 
set the number of threads for calls to BLAS to use

Input Parameter:
- `nt` - the number of threads

Options Database Key:
- `-blas_num_threads <nt>` - set the number of threads when PETSc is initialized

Level: intermediate

-seealso: `PetscInitialize()`, `PetscBLASGetNumThreads()`

# External Links
$(_doc_external("Sys/PetscBLASSetNumThreads"))
"""
function PetscBLASSetNumThreads(petsclib::PetscLibType, nt::PetscInt) end

@for_petsc function PetscBLASSetNumThreads(petsclib::$UnionPetscLib, nt::$PetscInt )

    @chk ccall(
               (:PetscBLASSetNumThreads, $petsc_library),
               PetscErrorCode,
               ($PetscInt,),
               nt,
              )


	return nothing
end 

"""
	nt::PetscInt = PetscBLASGetNumThreads(petsclib::PetscLibType) 
get the number of threads for calls to BLAS to use

Output Parameter:
- `nt` - the number of threads

Level: intermediate

-seealso: `PetscInitialize()`, `PetscBLASSetNumThreads()`

# External Links
$(_doc_external("Sys/PetscBLASGetNumThreads"))
"""
function PetscBLASGetNumThreads(petsclib::PetscLibType) end

@for_petsc function PetscBLASGetNumThreads(petsclib::$UnionPetscLib)
	nt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscBLASGetNumThreads, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscInt},),
               nt_,
              )

	nt = nt_[]

	return nt
end 

"""
	size::Csize_t = PetscDataTypeGetSize(petsclib::PetscLibType,ptype::PetscDataType) 
Gets the size (in bytes) of a PETSc datatype

Not Collective

Input Parameter:
- `ptype` - the PETSc datatype name (for example `PETSC_DOUBLE`)

Output Parameter:
- `size` - the size in bytes (for example the size of `PETSC_DOUBLE` is 8)

Level: advanced

-seealso: `PetscDataType`, `PetscDataTypeToMPIDataType()`

# External Links
$(_doc_external("Sys/PetscDataTypeGetSize"))
"""
function PetscDataTypeGetSize(petsclib::PetscLibType, ptype::PetscDataType) end

@for_petsc function PetscDataTypeGetSize(petsclib::$UnionPetscLib, ptype::PetscDataType )
	size_ = Ref{Csize_t}()

    @chk ccall(
               (:PetscDataTypeGetSize, $petsc_library),
               PetscErrorCode,
               (PetscDataType, Ptr{Csize_t}),
               ptype, size_,
              )

	size = size_[]

	return size
end 

"""
	ptype::PetscDataType,found::PetscBool = PetscDataTypeFromString(petsclib::PetscLibType,name::Vector{Cchar}) 
Gets the enum value of a PETSc datatype represented as a string

Not Collective

Input Parameter:
- `name` - the PETSc datatype name (for example, "double" or "real")

Output Parameters:
- `ptype` - the enum value, only valid if found is `PETSC_TRUE`
- `found` - the string matches one of the data types

Level: advanced

-seealso: `PetscDataType`, `PetscDataTypeToMPIDataType()`, `PetscDataTypeGetSize()`

# External Links
$(_doc_external("Sys/PetscDataTypeFromString"))
"""
function PetscDataTypeFromString(petsclib::PetscLibType, name::Vector{Cchar}) end

@for_petsc function PetscDataTypeFromString(petsclib::$UnionPetscLib, name::Vector{Cchar} )
	ptype_ = Ref{PetscDataType}()
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDataTypeFromString, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscDataType}, Ptr{PetscBool}),
               name, ptype_, found_,
              )

	ptype = unsafe_string(ptype_[])
	found = found_[]

	return ptype,found
end 

"""
	PetscKokkosInitializeCheck(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscKokkosInitializeCheck"))
"""
function PetscKokkosInitializeCheck(petsclib::PetscLibType) end

@for_petsc function PetscKokkosInitializeCheck(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscKokkosInitializeCheck, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscElementalInitializePackage(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscElementalInitializePackage"))
"""
function PetscElementalInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscElementalInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscElementalInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	isInitialized::PetscBool = PetscElementalInitialized(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscElementalInitialized"))
"""
function PetscElementalInitialized(petsclib::PetscLibType) end

@for_petsc function PetscElementalInitialized(petsclib::$UnionPetscLib)
	isInitialized_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscElementalInitialized, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               isInitialized_,
              )

	isInitialized = isInitialized_[]

	return isInitialized
end 

"""
	PetscElementalFinalizePackage(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscElementalFinalizePackage"))
"""
function PetscElementalFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscElementalFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscElementalFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	type::PetscMemType = PetscGetMemType(petsclib::PetscLibType,ptr::Cvoid) 
Query the `PetscMemType` of a pointer

Not Collective, No Fortran Support

Input Parameter:
- `ptr` - The pointer to query (may be `NULL`)

Output Parameter:
- `type` - The `PetscMemType` of the pointer

-seealso: `PetscMemType`, `PetscDeviceMalloc()`, `PetscDeviceCalloc()`, `PetscDeviceFree()`,
`PetscDeviceArrayCopy()`, `PetscDeviceArrayZero()`

# External Links
$(_doc_external("Sys/PetscGetMemType"))
"""
function PetscGetMemType(petsclib::PetscLibType, ptr::Cvoid) end

@for_petsc function PetscGetMemType(petsclib::$UnionPetscLib, ptr::Cvoid )
	type_ = Ref{PetscMemType}()

    @chk ccall(
               (:PetscGetMemType, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, Ptr{PetscMemType}),
               ptr, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscTimSort(petsclib::PetscLibType,n::PetscInt, arr::Cvoid, size::Csize_t, cmp::external, ctx::Cvoid) 
Sorts an array in place in increasing order using Tim Peters <https://bugs.python.org/file4451/timsort.txt> adaptive sorting algorithm.

Not Collective, No Fortran Support

Input Parameters:
- `n`    - number of values
- `arr`  - array to be sorted
- `size` - size in bytes of the datatype held in arr
- `cmp`  - function pointer to comparison function
- `ctx`  - optional context to be passed to comparison function, NULL if not needed

Output Parameter:
- `arr` - sorted array

Level: developer

-seealso: `PetscTimSortWithArray()`, `PetscIntSortSemiOrdered()`, `PetscRealSortSemiOrdered()`, `PetscMPIIntSortSemiOrdered()`

# External Links
$(_doc_external("Sys/PetscTimSort"))
"""
function PetscTimSort(petsclib::PetscLibType, n::PetscInt, arr::Cvoid, size::Csize_t, cmp::external, ctx::Cvoid) end

@for_petsc function PetscTimSort(petsclib::$UnionPetscLib, n::$PetscInt, arr::Cvoid, size::Csize_t, cmp::external, ctx::Cvoid )

    @chk ccall(
               (:PetscTimSort, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Cvoid}, Csize_t, external, Ptr{Cvoid}),
               n, arr, size, cmp, ctx,
              )


	return nothing
end 

"""
	PetscTimSortWithArray(petsclib::PetscLibType,n::PetscInt, arr::Cvoid, asize::Csize_t, barr::Cvoid, bsize::Csize_t, cmp::external, ctx::Cvoid) 
Sorts an array in place in increasing order using Tim Peters <https://bugs.python.org/file4451/timsort.txt> adaptive sorting algorithm and
reorders a second array to match the first. The arrays need not be the same type.

Not Collective, No Fortran Support

Input Parameters:
- `n`     - number of values
- `asize` - size in bytes of the datatype held in arr
- `bsize` - size in bytes of the datatype held in barr
- `cmp`   - function pointer to comparison function
- `ctx`   - optional context to be passed to comparison function, NULL if not needed

Input/Output Parameters:
- `arr`  - array to be sorted, on output it is sorted
- `barr` - array to be reordered, on output it is reordered

Level: developer

-seealso: `PetscTimSort()`, `PetscIntSortSemiOrderedWithArray()`, `PetscRealSortSemiOrderedWithArrayInt()`, `PetscMPIIntSortSemiOrderedWithArray()`

# External Links
$(_doc_external("Sys/PetscTimSortWithArray"))
"""
function PetscTimSortWithArray(petsclib::PetscLibType, n::PetscInt, arr::Cvoid, asize::Csize_t, barr::Cvoid, bsize::Csize_t, cmp::external, ctx::Cvoid) end

@for_petsc function PetscTimSortWithArray(petsclib::$UnionPetscLib, n::$PetscInt, arr::Cvoid, asize::Csize_t, barr::Cvoid, bsize::Csize_t, cmp::external, ctx::Cvoid )

    @chk ccall(
               (:PetscTimSortWithArray, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Cvoid}, Csize_t, Ptr{Cvoid}, Csize_t, external, Ptr{Cvoid}),
               n, arr, asize, barr, bsize, cmp, ctx,
              )


	return nothing
end 

"""
	PetscIntSortSemiOrdered(petsclib::PetscLibType,n::PetscInt, arr::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in increasing order.

Not Collective

Input Parameters:
- `n`   - number of values
- `arr` - array of integers

Output Parameter:
- `arr` - sorted array of integers

Level: intermediate

-seealso: `PetscTimSort()`, `PetscSortInt()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscIntSortSemiOrdered"))
"""
function PetscIntSortSemiOrdered(petsclib::PetscLibType, n::PetscInt, arr::Vector{PetscInt}) end

@for_petsc function PetscIntSortSemiOrdered(petsclib::$UnionPetscLib, n::$PetscInt, arr::Vector{$PetscInt} )

    @chk ccall(
               (:PetscIntSortSemiOrdered, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}),
               n, arr,
              )


	return nothing
end 

"""
	PetscIntSortSemiOrderedWithArray(petsclib::PetscLibType,n::PetscInt, arr1::Vector{PetscInt}, arr2::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in increasing order and reorders a second
`PetscInt` array to match the first.

Not Collective

Input Parameter:
- `n` - number of values

Input/Output Parameters:
- `arr1` - array of integers to be sorted, modified on output
- `arr2` - array of integers to be reordered, modified on output

Level: intermediate

-seealso: `PetscTimSortWithArray()`, `PetscSortIntWithArray()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscIntSortSemiOrderedWithArray"))
"""
function PetscIntSortSemiOrderedWithArray(petsclib::PetscLibType, n::PetscInt, arr1::Vector{PetscInt}, arr2::Vector{PetscInt}) end

@for_petsc function PetscIntSortSemiOrderedWithArray(petsclib::$UnionPetscLib, n::$PetscInt, arr1::Vector{$PetscInt}, arr2::Vector{$PetscInt} )

    @chk ccall(
               (:PetscIntSortSemiOrderedWithArray, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, arr1, arr2,
              )


	return nothing
end 

"""
	PetscMPIIntSortSemiOrdered(petsclib::PetscLibType,n::PetscInt, arr::Vector{PetscMPIInt}) 
Sorts an array of `PetscMPIInt` in place in increasing order.

Not Collective

Input Parameters:
- `n`   - number of values
- `arr` - array of `PetscMPIInt`

Output Parameter:
- `arr` - sorted array of integers

Level: intermediate

-seealso: `PetscTimSort()`, `PetscSortMPIInt()`

# External Links
$(_doc_external("Sys/PetscMPIIntSortSemiOrdered"))
"""
function PetscMPIIntSortSemiOrdered(petsclib::PetscLibType, n::PetscInt, arr::Vector{PetscMPIInt}) end

@for_petsc function PetscMPIIntSortSemiOrdered(petsclib::$UnionPetscLib, n::$PetscInt, arr::Vector{PetscMPIInt} )

    @chk ccall(
               (:PetscMPIIntSortSemiOrdered, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscMPIInt}),
               n, arr,
              )


	return nothing
end 

"""
	PetscMPIIntSortSemiOrderedWithArray(petsclib::PetscLibType,n::PetscInt, arr1::Vector{PetscMPIInt}, arr2::Vector{PetscMPIInt}) 
Sorts an array of `PetscMPIInt` in place in increasing order and reorders a second `PetscMPIInt`
array to match the first.

Not Collective

Input Parameter:
- `n` - number of values

Input/Output Parameters:
- `arr1` - array of integers to be sorted, modified on output
- `arr2` - array of integers to be reordered, modified on output

Level: intermediate

-seealso: `PetscTimSortWithArray()`, `PetscSortMPIIntWithArray()`, `PetscSortMPIIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscMPIIntSortSemiOrderedWithArray"))
"""
function PetscMPIIntSortSemiOrderedWithArray(petsclib::PetscLibType, n::PetscInt, arr1::Vector{PetscMPIInt}, arr2::Vector{PetscMPIInt}) end

@for_petsc function PetscMPIIntSortSemiOrderedWithArray(petsclib::$UnionPetscLib, n::$PetscInt, arr1::Vector{PetscMPIInt}, arr2::Vector{PetscMPIInt} )

    @chk ccall(
               (:PetscMPIIntSortSemiOrderedWithArray, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}),
               n, arr1, arr2,
              )


	return nothing
end 

"""
	PetscRealSortSemiOrdered(petsclib::PetscLibType,n::PetscInt, arr::Vector{PetscReal}) 
Sorts an array of `PetscReal` in place in increasing order.

Not Collective

Input Parameters:
- `n`   - number of values
- `arr` - array of `PetscReal`

Output Parameter:
- `arr` - sorted array of integers

Level: intermediate

-seealso: `PetscTimSort()`, `PetscSortReal()`, `PetscSortRealWithPermutation()`

# External Links
$(_doc_external("Sys/PetscRealSortSemiOrdered"))
"""
function PetscRealSortSemiOrdered(petsclib::PetscLibType, n::PetscInt, arr::Vector{PetscReal}) end

@for_petsc function PetscRealSortSemiOrdered(petsclib::$UnionPetscLib, n::$PetscInt, arr::Vector{$PetscReal} )

    @chk ccall(
               (:PetscRealSortSemiOrdered, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}),
               n, arr,
              )


	return nothing
end 

"""
	PetscRealSortSemiOrderedWithArrayInt(petsclib::PetscLibType,n::PetscInt, arr1::Vector{PetscReal}, arr2::Vector{PetscInt}) 
Sorts an array of `PetscReal` in place in increasing order and reorders a second
array of `PetscInt` to match the first.

Not Collective

Input Parameter:
- `n` - number of values

Input/Output Parameters:
- `arr1` - array of `PetscReal` to be sorted, modified on output
- `arr2` - array of `PetscInt` to be reordered, modified on output

Level: intermediate

-seealso: `PetscTimSortWithArray()`, `PetscSortRealWithArrayInt()`, `PetscSortRealWithPermutation()`

# External Links
$(_doc_external("Sys/PetscRealSortSemiOrderedWithArrayInt"))
"""
function PetscRealSortSemiOrderedWithArrayInt(petsclib::PetscLibType, n::PetscInt, arr1::Vector{PetscReal}, arr2::Vector{PetscInt}) end

@for_petsc function PetscRealSortSemiOrderedWithArrayInt(petsclib::$UnionPetscLib, n::$PetscInt, arr1::Vector{$PetscReal}, arr2::Vector{$PetscInt} )

    @chk ccall(
               (:PetscRealSortSemiOrderedWithArrayInt, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
               n, arr1, arr2,
              )


	return nothing
end 

"""
	slope::PetscReal,intercept::PetscReal = PetscLinearRegression(petsclib::PetscLibType,n::PetscInt, x::Vector{PetscReal}, y::Vector{PetscReal}) 
Gives the best least

Input Parameters:
- `n` - The number of points
- `x` - The x-values
- `y` - The y-values

Output Parameters:
- `slope`     - The slope of the best-fit line
- `intercept` - The y-intercept of the best-fit line

Level: intermediate

-seealso: `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Sys/PetscLinearRegression"))
"""
function PetscLinearRegression(petsclib::PetscLibType, n::PetscInt, x::Vector{PetscReal}, y::Vector{PetscReal}) end

@for_petsc function PetscLinearRegression(petsclib::$UnionPetscLib, n::$PetscInt, x::Vector{$PetscReal}, y::Vector{$PetscReal} )
	slope_ = Ref{$PetscReal}()
	intercept_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscLinearRegression, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               n, x, y, slope_, intercept_,
              )

	slope = slope_[]
	intercept = intercept_[]

	return slope,intercept
end 

"""
	sorted::PetscBool = PetscSortedInt(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}) 
Determines whether the `PetscInt` array is sorted.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`

Output Parameter:
- `sorted` - flag whether the array is sorted

Level: intermediate

-seealso: `PetscSortInt()`, `PetscSortedMPIInt()`, `PetscSortedReal()`

# External Links
$(_doc_external("Sys/PetscSortedInt"))
"""
function PetscSortedInt(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}) end

@for_petsc function PetscSortedInt(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt} )
	sorted_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSortedInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, X, sorted_,
              )

	sorted = sorted_[]

	return sorted
end 

"""
	sorted::PetscBool = PetscSortedInt64(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt64}) 
Determines whether the `PetscInt64` array is sorted.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt64`

Output Parameter:
- `sorted` - flag whether the array is sorted

Level: intermediate

-seealso: `PetscSortInt64()`, `PetscSortInt()`, `PetscSortedMPIInt()`, `PetscSortedReal()`

# External Links
$(_doc_external("Sys/PetscSortedInt64"))
"""
function PetscSortedInt64(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt64}) end

@for_petsc function PetscSortedInt64(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt64} )
	sorted_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSortedInt64, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt64}, Ptr{PetscBool}),
               n, X, sorted_,
              )

	sorted = sorted_[]

	return sorted
end 

"""
	PetscSortInt(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in increasing order.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`

-seealso: `PetscIntSortSemiOrdered()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortInt"))
"""
function PetscSortInt(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}) end

@for_petsc function PetscSortInt(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}),
               n, X,
              )


	return nothing
end 

"""
	PetscSortInt64(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt64}) 
Sorts an array of `PetscInt64` in place in increasing order.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt64`

-seealso: `PetscSortInt()`

# External Links
$(_doc_external("Sys/PetscSortInt64"))
"""
function PetscSortInt64(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt64}) end

@for_petsc function PetscSortInt64(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt64} )

    @chk ccall(
               (:PetscSortInt64, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt64}),
               n, X,
              )


	return nothing
end 

"""
	PetscSortCount(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscCount}) 
Sorts an array of `PetscCount` in place in increasing order.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscCount`

-seealso: `PetscSortInt()`

# External Links
$(_doc_external("Sys/PetscSortCount"))
"""
function PetscSortCount(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscCount}) end

@for_petsc function PetscSortCount(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{PetscCount} )

    @chk ccall(
               (:PetscSortCount, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscCount}),
               n, X,
              )


	return nothing
end 

"""
	PetscSortReverseInt(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in decreasing order.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`

Level: intermediate

-seealso: `PetscIntSortSemiOrdered()`, `PetscSortInt()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortReverseInt"))
"""
function PetscSortReverseInt(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}) end

@for_petsc function PetscSortReverseInt(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortReverseInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}),
               n, X,
              )


	return nothing
end 

"""
	PetscSortedRemoveDupsInt(petsclib::PetscLibType,n::PetscInt, X::Vector{PetscInt}) 
Removes all duplicate entries of a sorted `PetscInt` array

Not Collective

Input Parameters:
- `n` - number of values
- `X` - sorted array of `PetscInt`

Output Parameter:
- `n` - number of non-redundant values

Level: intermediate

-seealso: `PetscSortInt()`

# External Links
$(_doc_external("Sys/PetscSortedRemoveDupsInt"))
"""
function PetscSortedRemoveDupsInt(petsclib::PetscLibType, n::PetscInt, X::Vector{PetscInt}) end

@for_petsc function PetscSortedRemoveDupsInt(petsclib::$UnionPetscLib, n::$PetscInt, X::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortedRemoveDupsInt, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, X,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscSortedCheckDupsInt(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}) 
Checks if a sorted `PetscInt` array has duplicates

Not Collective

Input Parameters:
- `n` - number of values
- `X` - sorted array of `PetscInt`

Output Parameter:
- `flg` - True if the array has duplications, otherwise false

Level: intermediate

-seealso: `PetscSortInt()`, `PetscCheckDupsInt()`, `PetscSortRemoveDupsInt()`, `PetscSortedRemoveDupsInt()`

# External Links
$(_doc_external("Sys/PetscSortedCheckDupsInt"))
"""
function PetscSortedCheckDupsInt(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}) end

@for_petsc function PetscSortedCheckDupsInt(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSortedCheckDupsInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, X, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PetscSortedCheckDupsCount(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscCount}) 
Checks if a sorted `PetscCount` array has duplicates

Not Collective

Input Parameters:
- `n` - number of values
- `X` - sorted array of `PetscCount`

Output Parameter:
- `flg` - True if the array has duplications, otherwise false

Level: intermediate

-seealso: `PetscCount`, `PetscSortCount()`, `PetscSortedCheckDupsInt()`

# External Links
$(_doc_external("Sys/PetscSortedCheckDupsCount"))
"""
function PetscSortedCheckDupsCount(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscCount}) end

@for_petsc function PetscSortedCheckDupsCount(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{PetscCount} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSortedCheckDupsCount, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscCount}, Ptr{PetscBool}),
               n, X, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscSortRemoveDupsInt(petsclib::PetscLibType,n::PetscInt, X::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in increasing order removes all duplicate entries

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`

Output Parameter:
- `n` - number of non-redundant values

Level: intermediate

-seealso: `PetscIntSortSemiOrdered()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortedRemoveDupsInt()`

# External Links
$(_doc_external("Sys/PetscSortRemoveDupsInt"))
"""
function PetscSortRemoveDupsInt(petsclib::PetscLibType, n::PetscInt, X::Vector{PetscInt}) end

@for_petsc function PetscSortRemoveDupsInt(petsclib::$UnionPetscLib, n::$PetscInt, X::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortRemoveDupsInt, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, X,
              )


	return nothing
end 

"""
	loc::PetscInt = PetscFindInt(petsclib::PetscLibType,key::PetscInt, n::PetscCount, X::Vector{PetscInt}) 
Finds the location of a `PetscInt` key in a sorted array of `PetscInt`

Not Collective

Input Parameters:
- `key` - the `PetscInt` key to locate
- `n`   - number of values in the array
- `X`   - array of `PetscInt`

Output Parameter:
- `loc` - the location if found, otherwise -(slot+1) where slot is the place the value would go

Level: intermediate

-seealso: `PetscIntSortSemiOrdered()`, `PetscSortInt()`, `PetscSortIntWithArray()`, `PetscSortRemoveDupsInt()`

# External Links
$(_doc_external("Sys/PetscFindInt"))
"""
function PetscFindInt(petsclib::PetscLibType, key::PetscInt, n::PetscCount, X::Vector{PetscInt}) end

@for_petsc function PetscFindInt(petsclib::$UnionPetscLib, key::$PetscInt, n::PetscCount, X::Vector{$PetscInt} )
	loc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFindInt, $petsc_library),
               PetscErrorCode,
               ($PetscInt, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
               key, n, X, loc_,
              )

	loc = loc_[]

	return loc
end 

"""
	PetscFindCount(petsclib::PetscLibType,key::PetscCount, n::PetscCount, X::Vector{PetscCount}, loc::PetscCount) 
Finds the location of a `PetscCount` key in a sorted array of `PetscCount`

Not Collective

Input Parameters:
- `key` - the `PetscCount` key to locate
- `n`   - number of values in the array
- `X`   - array of `PetscCount`

Output Parameter:
- `loc` - the location if found, otherwise -(slot+1) where slot is the place the value would go

Level: intermediate

-seealso: `PetscCount`, `PetscSortCount()`

# External Links
$(_doc_external("Sys/PetscFindCount"))
"""
function PetscFindCount(petsclib::PetscLibType, key::PetscCount, n::PetscCount, X::Vector{PetscCount}, loc::PetscCount) end

@for_petsc function PetscFindCount(petsclib::$UnionPetscLib, key::PetscCount, n::PetscCount, X::Vector{PetscCount}, loc::PetscCount )

    @chk ccall(
               (:PetscFindCount, $petsc_library),
               PetscErrorCode,
               (PetscCount, PetscCount, Ptr{PetscCount}, Ptr{PetscCount}),
               key, n, X, loc,
              )


	return nothing
end 

"""
	dups::PetscBool = PetscCheckDupsInt(petsclib::PetscLibType,n::PetscInt, X::Vector{PetscInt}) 
Checks if an `PetscInt` array has duplicates

Not Collective

Input Parameters:
- `n` - number of values in the array
- `X` - array of `PetscInt`

Output Parameter:
- `dups` - True if the array has dups, otherwise false

Level: intermediate

-seealso: `PetscSortRemoveDupsInt()`, `PetscSortedCheckDupsInt()`

# External Links
$(_doc_external("Sys/PetscCheckDupsInt"))
"""
function PetscCheckDupsInt(petsclib::PetscLibType, n::PetscInt, X::Vector{PetscInt}) end

@for_petsc function PetscCheckDupsInt(petsclib::$UnionPetscLib, n::$PetscInt, X::Vector{$PetscInt} )
	dups_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscCheckDupsInt, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, X, dups_,
              )

	dups = dups_[]

	return dups
end 

"""
	loc::PetscInt = PetscFindMPIInt(petsclib::PetscLibType,key::PetscMPIInt, n::PetscCount, X::Vector{PetscMPIInt}) 
Finds `PetscMPIInt` in a sorted array of `PetscMPIInt`

Not Collective

Input Parameters:
- `key` - the integer to locate
- `n`   - number of values in the array
- `X`   - array of `PetscMPIInt`

Output Parameter:
- `loc` - the location if found, otherwise -(slot+1) where slot is the place the value would go

Level: intermediate

-seealso: `PetscMPIIntSortSemiOrdered()`, `PetscSortInt()`, `PetscSortIntWithArray()`, `PetscSortRemoveDupsInt()`

# External Links
$(_doc_external("Sys/PetscFindMPIInt"))
"""
function PetscFindMPIInt(petsclib::PetscLibType, key::PetscMPIInt, n::PetscCount, X::Vector{PetscMPIInt}) end

@for_petsc function PetscFindMPIInt(petsclib::$UnionPetscLib, key::PetscMPIInt, n::PetscCount, X::Vector{PetscMPIInt} )
	loc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFindMPIInt, $petsc_library),
               PetscErrorCode,
               (PetscMPIInt, PetscCount, Ptr{PetscMPIInt}, Ptr{$PetscInt}),
               key, n, X, loc_,
              )

	loc = loc_[]

	return loc
end 

"""
	PetscSortIntWithArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in increasing order;
changes a second array of `PetscInt` to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`
- `Y` - second array of `PetscInt`

Level: intermediate

-seealso: `PetscIntSortSemiOrderedWithArray()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithCountArray()`

# External Links
$(_doc_external("Sys/PetscSortIntWithArray"))
"""
function PetscSortIntWithArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscInt}) end

@for_petsc function PetscSortIntWithArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortIntWithArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, X, Y,
              )


	return nothing
end 

"""
	PetscSortIntWithArrayPair(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscInt}, Z::Vector{PetscInt}) 
Sorts an array of `PetscInt` in place in increasing order;
changes a pair of `PetscInt` arrays to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PestcInt`
- `Y` - second array of `PestcInt` (first array of the pair)
- `Z` - third array of `PestcInt` (second array of the pair)

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortIntWithArray()`, `PetscIntSortSemiOrdered()`, `PetscSortIntWithIntCountArrayPair()`

# External Links
$(_doc_external("Sys/PetscSortIntWithArrayPair"))
"""
function PetscSortIntWithArrayPair(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscInt}, Z::Vector{PetscInt}) end

@for_petsc function PetscSortIntWithArrayPair(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Vector{$PetscInt}, Z::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortIntWithArrayPair, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, X, Y, Z,
              )


	return nothing
end 

"""
	PetscSortIntWithMPIIntArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscMPIInt}) 
Sorts an array of `PetscInt` in place in increasing order;
changes a second array of `PetscMPI` to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`
- `Y` - second array of `PetscMPIInt`

Level: intermediate

-seealso: `PetscIntSortSemiOrderedWithArray()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscSortIntWithMPIIntArray"))
"""
function PetscSortIntWithMPIIntArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscMPIInt}) end

@for_petsc function PetscSortIntWithMPIIntArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Vector{PetscMPIInt} )

    @chk ccall(
               (:PetscSortIntWithMPIIntArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{PetscMPIInt}),
               n, X, Y,
              )


	return nothing
end 

"""
	PetscSortIntWithCountArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscCount}) 
Sorts an array of `PetscInt` in place in increasing order;
changes a second array of `PetscCount` to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`
- `Y` - second array of `PetscCount`

Level: intermediate

-seealso: `PetscIntSortSemiOrderedWithArray()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscSortIntWithCountArray"))
"""
function PetscSortIntWithCountArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscCount}) end

@for_petsc function PetscSortIntWithCountArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Vector{PetscCount} )

    @chk ccall(
               (:PetscSortIntWithCountArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{PetscCount}),
               n, X, Y,
              )


	return nothing
end 

"""
	PetscSortIntWithIntCountArrayPair(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscInt}, Z::Vector{PetscCount}) 
Sorts an array of `PetscInt` in place in increasing order;
changes a `PetscInt`  array and a `PetscCount` array to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`
- `Y` - second array of `PetscInt` (first array of the pair)
- `Z` - third array of `PetscCount` (second array of the pair)

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortIntWithArray()`, `PetscIntSortSemiOrdered()`, `PetscSortIntWithArrayPair()`

# External Links
$(_doc_external("Sys/PetscSortIntWithIntCountArrayPair"))
"""
function PetscSortIntWithIntCountArrayPair(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscInt}, Z::Vector{PetscCount}) end

@for_petsc function PetscSortIntWithIntCountArrayPair(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Vector{$PetscInt}, Z::Vector{PetscCount} )

    @chk ccall(
               (:PetscSortIntWithIntCountArrayPair, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{PetscCount}),
               n, X, Y, Z,
              )


	return nothing
end 

"""
	sorted::PetscBool = PetscSortedMPIInt(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscMPIInt}) 
Determines whether the `PetscMPIInt` array is sorted.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscMPIInt`

Output Parameter:
- `sorted` - flag whether the array is sorted

Level: intermediate

-seealso: `PetscMPIIntSortSemiOrdered()`, `PetscSortMPIInt()`, `PetscSortedInt()`, `PetscSortedReal()`

# External Links
$(_doc_external("Sys/PetscSortedMPIInt"))
"""
function PetscSortedMPIInt(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscMPIInt}) end

@for_petsc function PetscSortedMPIInt(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{PetscMPIInt} )
	sorted_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSortedMPIInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscMPIInt}, Ptr{PetscBool}),
               n, X, sorted_,
              )

	sorted = sorted_[]

	return sorted
end 

"""
	PetscSortMPIInt(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscMPIInt}) 
Sorts an array of `PetscMPIInt` in place in increasing order.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscMPIInt`

Level: intermediate

-seealso: `PetscMPIIntSortSemiOrdered()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortMPIInt"))
"""
function PetscSortMPIInt(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscMPIInt}) end

@for_petsc function PetscSortMPIInt(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{PetscMPIInt} )

    @chk ccall(
               (:PetscSortMPIInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscMPIInt}),
               n, X,
              )


	return nothing
end 

"""
	PetscSortRemoveDupsMPIInt(petsclib::PetscLibType,n::PetscInt, X::Vector{PetscMPIInt}) 
Sorts an array of `PetscMPIInt` in place in increasing order removes all duplicate entries

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscMPIInt`

Output Parameter:
- `n` - number of non-redundant values

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`

# External Links
$(_doc_external("Sys/PetscSortRemoveDupsMPIInt"))
"""
function PetscSortRemoveDupsMPIInt(petsclib::PetscLibType, n::PetscInt, X::Vector{PetscMPIInt}) end

@for_petsc function PetscSortRemoveDupsMPIInt(petsclib::$UnionPetscLib, n::$PetscInt, X::Vector{PetscMPIInt} )

    @chk ccall(
               (:PetscSortRemoveDupsMPIInt, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscInt}, Ptr{PetscMPIInt}),
               n, X,
              )


	return nothing
end 

"""
	PetscSortMPIIntWithArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscMPIInt}, Y::Vector{PetscMPIInt}) 
Sorts an array of `PetscMPIInt` in place in increasing order;
changes a second `PetscMPIInt` array to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscMPIInt`
- `Y` - second array of `PetscMPIInt`

Level: intermediate

-seealso: `PetscMPIIntSortSemiOrderedWithArray()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`

# External Links
$(_doc_external("Sys/PetscSortMPIIntWithArray"))
"""
function PetscSortMPIIntWithArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscMPIInt}, Y::Vector{PetscMPIInt}) end

@for_petsc function PetscSortMPIIntWithArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{PetscMPIInt}, Y::Vector{PetscMPIInt} )

    @chk ccall(
               (:PetscSortMPIIntWithArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}),
               n, X, Y,
              )


	return nothing
end 

"""
	PetscSortMPIIntWithIntArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscMPIInt}, Y::Vector{PetscInt}) 
Sorts an array of `PetscMPIInt` in place in increasing order;
changes a second array of `PetscInt` to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscMPIInt`
- `Y` - second array of `PetscInt`

Level: intermediate

-seealso: `PetscSortMPIIntWithArray()`, `PetscIntSortSemiOrderedWithArray()`, `PetscTimSortWithArray()`

# External Links
$(_doc_external("Sys/PetscSortMPIIntWithIntArray"))
"""
function PetscSortMPIIntWithIntArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscMPIInt}, Y::Vector{PetscInt}) end

@for_petsc function PetscSortMPIIntWithIntArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{PetscMPIInt}, Y::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortMPIIntWithIntArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscMPIInt}, Ptr{$PetscInt}),
               n, X, Y,
              )


	return nothing
end 

"""
	PetscSortIntWithScalarArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscScalar}) 
Sorts an array of `PetscInt` in place in increasing order;
changes a second `PetscScalar` array to match the sorted first array.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of `PetscInt`
- `Y` - second array of `PetscScalar`

Level: intermediate

-seealso: `PetscTimSortWithArray()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscSortIntWithScalarArray"))
"""
function PetscSortIntWithScalarArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Vector{PetscScalar}) end

@for_petsc function PetscSortIntWithScalarArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscSortIntWithScalarArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               n, X, Y,
              )


	return nothing
end 

"""
	PetscSortIntWithDataArray(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscInt}, Y::Cvoid, size::Csize_t, t2::Cvoid) 
Sorts an array of `PetscInt` in place in increasing order;
changes a second array to match the sorted first INTEGER array.  Unlike other sort routines, the user must
provide workspace (the size of an element in the data array) to use when sorting.

Not Collective, No Fortran Support

Input Parameters:
- `n`    - number of values
- `X`    - array of `PetscInt`
- `Y`    - second array of data
- `size` - sizeof elements in the data array in bytes
- `t2`   - workspace of "size" bytes used when sorting

Level: intermediate

-seealso: `PetscTimSortWithArray()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscSortIntWithDataArray"))
"""
function PetscSortIntWithDataArray(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscInt}, Y::Cvoid, size::Csize_t, t2::Cvoid) end

@for_petsc function PetscSortIntWithDataArray(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscInt}, Y::Cvoid, size::Csize_t, t2::Cvoid )

    @chk ccall(
               (:PetscSortIntWithDataArray, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscInt}, Ptr{Cvoid}, Csize_t, Ptr{Cvoid}),
               n, X, Y, size, t2,
              )


	return nothing
end 

"""
	n::PetscInt,L::PetscInt = PetscMergeIntArray(petsclib::PetscLibType,an::PetscInt, aI::Vector{PetscInt}, bn::PetscInt, bI::Vector{PetscInt}) 
Merges two SORTED `PetscInt` arrays, removes duplicate elements.

Not Collective

Input Parameters:
- `an` - number of values in the first array
- `aI` - first sorted array of `PetscInt`
- `bn` - number of values in the second array
- `bI` - second array of `PetscInt`

Output Parameters:
- `n` - number of values in the merged array
- `L` - merged sorted array, this is allocated if an array is not provided

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscMergeIntArray"))
"""
function PetscMergeIntArray(petsclib::PetscLibType, an::PetscInt, aI::Vector{PetscInt}, bn::PetscInt, bI::Vector{PetscInt}) end

@for_petsc function PetscMergeIntArray(petsclib::$UnionPetscLib, an::$PetscInt, aI::Vector{$PetscInt}, bn::$PetscInt, bI::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}()
	L_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscMergeIntArray, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt),
               an, aI, bn, bI, n_, L_,
              )

	n = n_[]
	L = L_[]

	return n,L
end 

"""
	n::PetscInt,L::Vector{PetscInt},J::Vector{PetscInt} = PetscMergeIntArrayPair(petsclib::PetscLibType,an::PetscInt, aI::Vector{PetscInt}, aJ::Vector{PetscInt}, bn::PetscInt, bI::Vector{PetscInt}, bJ::Vector{PetscInt}) 
Merges two SORTED `PetscInt` arrays that share NO common values along with an additional array of `PetscInt`.
The additional arrays are the same length as sorted arrays and are merged
in the order determined by the merging of the sorted pair.

Not Collective

Input Parameters:
- `an` - number of values in the first array
- `aI` - first sorted array of `PetscInt`
- `aJ` - first additional array of `PetscInt`
- `bn` - number of values in the second array
- `bI` - second array of `PetscInt`
- `bJ` - second additional of `PetscInt`

Output Parameters:
- `n` - number of values in the merged array (== an + bn)
- `L` - merged sorted array
- `J` - merged additional array

-seealso: `PetscIntSortSemiOrdered()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscMergeIntArrayPair"))
"""
function PetscMergeIntArrayPair(petsclib::PetscLibType, an::PetscInt, aI::Vector{PetscInt}, aJ::Vector{PetscInt}, bn::PetscInt, bI::Vector{PetscInt}, bJ::Vector{PetscInt}) end

@for_petsc function PetscMergeIntArrayPair(petsclib::$UnionPetscLib, an::$PetscInt, aI::Vector{$PetscInt}, aJ::Vector{$PetscInt}, bn::$PetscInt, bI::Vector{$PetscInt}, bJ::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}()
	L_ = Ref{Ptr{$PetscInt}}()
	J_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscMergeIntArrayPair, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               an, aI, aJ, bn, bI, bJ, n_, L_, J_,
              )

	n = n_[]
	L = unsafe_wrap(Array, L_[], VecGetLocalSize(petsclib, x); own = false)
	J = unsafe_wrap(Array, J_[], VecGetLocalSize(petsclib, x); own = false)

	return n,L,J
end 

"""
	n::PetscInt = PetscMergeMPIIntArray(petsclib::PetscLibType,an::PetscInt, aI::Vector{PetscMPIInt}, bn::PetscInt, bI::Vector{PetscMPIInt}, L::PetscMPIInt) 
Merges two SORTED `PetscMPIInt` arrays.

Not Collective

Input Parameters:
- `an` - number of values in the first array
- `aI` - first sorted array of `PetscMPIInt`
- `bn` - number of values in the second array
- `bI` - second array of `PetscMPIInt`

Output Parameters:
- `n` - number of values in the merged array (<= an + bn)
- `L` - merged sorted array, allocated if address of NULL pointer is passed

Level: intermediate

-seealso: `PetscIntSortSemiOrdered()`, `PetscSortReal()`, `PetscSortIntWithPermutation()`, `PetscSortInt()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscMergeMPIIntArray"))
"""
function PetscMergeMPIIntArray(petsclib::PetscLibType, an::PetscInt, aI::Vector{PetscMPIInt}, bn::PetscInt, bI::Vector{PetscMPIInt}, L::PetscMPIInt) end

@for_petsc function PetscMergeMPIIntArray(petsclib::$UnionPetscLib, an::$PetscInt, aI::Vector{PetscMPIInt}, bn::$PetscInt, bI::Vector{PetscMPIInt}, L::PetscMPIInt )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscMergeMPIIntArray, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscMPIInt}, $PetscInt, Ptr{PetscMPIInt}, Ptr{$PetscInt}, PetscMPIInt),
               an, aI, bn, bI, n_, L,
              )

	n = n_[]

	return n
end 

"""
	Nlevels::PetscInt,Level::PetscInt,Levelcnt::PetscInt,Idbylevel::PetscInt,Column::PetscInt = PetscProcessTree(petsclib::PetscLibType,n::PetscInt, mask::Vector{PetscBool}, parentid::Vector{PetscInt}) 
Prepares tree data to be displayed graphically

Not Collective, No Fortran Support

Input Parameters:
- `n`        - number of values
- `mask`     - indicates those entries in the tree, location 0 is always masked
- `parentid` - indicates the parent of each entry

Output Parameters:
- `Nlevels`   - the number of levels
- `Level`     - for each node tells its level
- `Levelcnt`  - the number of nodes on each level
- `Idbylevel` - a list of ids on each of the levels, first level followed by second etc
- `Column`    - for each id tells its column index

Level: developer

-seealso: `PetscSortReal()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscProcessTree"))
"""
function PetscProcessTree(petsclib::PetscLibType, n::PetscInt, mask::Vector{PetscBool}, parentid::Vector{PetscInt}) end

@for_petsc function PetscProcessTree(petsclib::$UnionPetscLib, n::$PetscInt, mask::Vector{PetscBool}, parentid::Vector{$PetscInt} )
	Nlevels_ = Ref{$PetscInt}()
	Level_ = Ref{$PetscInt}()
	Levelcnt_ = Ref{$PetscInt}()
	Idbylevel_ = Ref{$PetscInt}()
	Column_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscProcessTree, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscBool}, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
               n, mask, parentid, Nlevels_, Level_, Levelcnt_, Idbylevel_, Column_,
              )

	Nlevels = Nlevels_[]
	Level = Level_[]
	Levelcnt = Levelcnt_[]
	Idbylevel = Idbylevel_[]
	Column = Column_[]

	return Nlevels,Level,Levelcnt,Idbylevel,Column
end 

"""
	is_sorted::PetscBool = PetscParallelSortedInt(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, keys::Vector{PetscInt}) 
Check whether a `PetscInt` array, distributed over a communicator, is globally sorted.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `n`    - the local number of `PetscInt`
- `keys` - the local array of `PetscInt`

Output Parameters:
- `is_sorted` - whether the array is globally sorted

Level: developer

-seealso: `PetscParallelSortInt()`

# External Links
$(_doc_external("Sys/PetscParallelSortedInt"))
"""
function PetscParallelSortedInt(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, keys::Vector{PetscInt}) end

@for_petsc function PetscParallelSortedInt(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, keys::Vector{$PetscInt} )
	is_sorted_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscParallelSortedInt, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               comm, n, keys, is_sorted_,
              )

	is_sorted = is_sorted_[]

	return is_sorted
end 

"""
	PetscCommBuildTwoSidedSetType(petsclib::PetscLibType,comm::MPI_Comm, twosided::PetscBuildTwoSidedType) 
set algorithm to use when building two

Logically Collective

Input Parameters:
- `comm`     - `PETSC_COMM_WORLD`
- `twosided` - algorithm to use in subsequent calls to `PetscCommBuildTwoSided()`

Level: developer

-seealso: `PetscCommBuildTwoSided()`, `PetscCommBuildTwoSidedGetType()`, `PetscBuildTwoSidedType`

# External Links
$(_doc_external("Sys/PetscCommBuildTwoSidedSetType"))
"""
function PetscCommBuildTwoSidedSetType(petsclib::PetscLibType, comm::MPI_Comm, twosided::PetscBuildTwoSidedType) end

@for_petsc function PetscCommBuildTwoSidedSetType(petsclib::$UnionPetscLib, comm::MPI_Comm, twosided::PetscBuildTwoSidedType )

    @chk ccall(
               (:PetscCommBuildTwoSidedSetType, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscBuildTwoSidedType),
               comm, twosided,
              )


	return nothing
end 

"""
	comm::MPI_Comm,twosided::PetscBuildTwoSidedType = PetscCommBuildTwoSidedGetType(petsclib::PetscLibType) 
get algorithm used when building two

Logically Collective

Output Parameters:
- `comm`     - communicator on which to query algorithm
- `twosided` - algorithm to use for `PetscCommBuildTwoSided()`

Level: developer

-seealso: `PetscCommBuildTwoSided()`, `PetscCommBuildTwoSidedSetType()`, `PetscBuildTwoSidedType`

# External Links
$(_doc_external("Sys/PetscCommBuildTwoSidedGetType"))
"""
function PetscCommBuildTwoSidedGetType(petsclib::PetscLibType) end

@for_petsc function PetscCommBuildTwoSidedGetType(petsclib::$UnionPetscLib)
	comm_ = Ref{MPI_Comm}()
	twosided_ = Ref{PetscBuildTwoSidedType}()

    @chk ccall(
               (:PetscCommBuildTwoSidedGetType, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscBuildTwoSidedType}),
               comm_, twosided_,
              )

	comm = comm_[]
	twosided = unsafe_string(twosided_[])

	return comm,twosided
end 

"""
	PetscSplitOwnershipBlock(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt) 
Given a global (or local) length determines a local
(or global) length via a simple formula. Splits so each processors local size
is divisible by the block size.

Collective (if `N` is `PETSC_DECIDE`)

Input Parameters:
- `comm` - MPI communicator that shares the object being divided
- `bs`   - block size
- `n`    - local length (or `PETSC_DECIDE` to have it set)
- `N`    - global length (or `PETSC_DECIDE`)

Level: developer

-seealso: `PetscSplitOwnership()`, `PetscSplitOwnershipEqual()`

# External Links
$(_doc_external("Sys/PetscSplitOwnershipBlock"))
"""
function PetscSplitOwnershipBlock(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt) end

@for_petsc function PetscSplitOwnershipBlock(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt )

    @chk ccall(
               (:PetscSplitOwnershipBlock, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               comm, bs, n, N,
              )


	return nothing
end 

"""
	PetscSplitOwnership(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Given a global (or local) length determines a local
(or global) length via a simple formula

Collective (if `n` or `N` is `PETSC_DECIDE` or `PETSC_DETERMINE`)

Input Parameters:
- `comm` - MPI communicator that shares the object being divided
- `n`    - local length (or `PETSC_DECIDE` to have it set)
- `N`    - global length (or `PETSC_DETERMINE` to have it set)

Level: developer

-seealso: `PetscSplitOwnershipBlock()`, `PetscSplitOwnershipEqual()`, `PETSC_DECIDE`, `PETSC_DETERMINE`

# External Links
$(_doc_external("Sys/PetscSplitOwnership"))
"""
function PetscSplitOwnership(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function PetscSplitOwnership(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )

    @chk ccall(
               (:PetscSplitOwnership, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}),
               comm, n, N,
              )


	return nothing
end 

"""
	PetscSplitOwnershipEqual(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Given a global (or local) length determines a local
(or global) length via a simple formula, trying to have all local lengths equal

Collective (if `n` or `N` is `PETSC_DECIDE`)

Input Parameters:
- `comm` - MPI communicator that shares the object being divided
- `n`    - local length (or `PETSC_DECIDE` to have it set)
- `N`    - global length (or `PETSC_DECIDE`)

Level: developer

-seealso: `PetscSplitOwnership()`, `PetscSplitOwnershipBlock()`

# External Links
$(_doc_external("Sys/PetscSplitOwnershipEqual"))
"""
function PetscSplitOwnershipEqual(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function PetscSplitOwnershipEqual(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )

    @chk ccall(
               (:PetscSplitOwnershipEqual, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}),
               comm, n, N,
              )


	return nothing
end 

"""
	PetscStrToArray(petsclib::PetscLibType,s::Vector{Cchar}, sp::Cchar, argc::Cint, args::Cchar) 
Separates a string by a character (for example ' ' or '\n') and creates an array of strings

Not Collective; No Fortran Support

Input Parameters:
- `s`  - pointer to string
- `sp` - separator character

Output Parameters:
- `argc` - the number of entries in `args`
- `args` - an array of the entries with a `NULL` at the end

Level: intermediate

-seealso: `PetscStrToArrayDestroy()`, `PetscToken`, `PetscTokenCreate()`

# External Links
$(_doc_external("Sys/PetscStrToArray"))
"""
function PetscStrToArray(petsclib::PetscLibType, s::Vector{Cchar}, sp::Cchar, argc::Cint, args::Cchar) end

@for_petsc function PetscStrToArray(petsclib::$UnionPetscLib, s::Vector{Cchar}, sp::Cchar, argc::Cint, args::Cchar )

    @chk ccall(
               (:PetscStrToArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{Cint}, Cchar),
               s, sp, argc, args,
              )


	return nothing
end 

"""
	PetscStrToArrayDestroy(petsclib::PetscLibType,argc::Cint, args::Cchar) 
Frees array created with `PetscStrToArray()`.

Not Collective; No Fortran Support

Output Parameters:
- `argc` - the number of arguments
- `args` - the array of arguments

Level: intermediate

-seealso: `PetscStrToArray()`

# External Links
$(_doc_external("Sys/PetscStrToArrayDestroy"))
"""
function PetscStrToArrayDestroy(petsclib::PetscLibType, argc::Cint, args::Cchar) end

@for_petsc function PetscStrToArrayDestroy(petsclib::$UnionPetscLib, argc::Cint, args::Cchar )

    @chk ccall(
               (:PetscStrToArrayDestroy, $petsc_library),
               PetscErrorCode,
               (Cint, Cchar),
               argc, args,
              )


	return nothing
end 

"""
	PetscStrArrayallocpy(petsclib::PetscLibType,list::Cchar, t::Cchar) 
Allocates space to hold a copy of an array of strings then copies the strings

Not Collective; No Fortran Support

Input Parameter:
- `list` - pointer to array of strings (final string is a `NULL`)

Output Parameter:
- `t` - the copied array string

Level: intermediate

-seealso: `PetscStrallocpy()`, `PetscStrArrayDestroy()`, `PetscStrNArrayallocpy()`

# External Links
$(_doc_external("Sys/PetscStrArrayallocpy"))
"""
function PetscStrArrayallocpy(petsclib::PetscLibType, list::Cchar, t::Cchar) end

@for_petsc function PetscStrArrayallocpy(petsclib::$UnionPetscLib, list::Cchar, t::Cchar )

    @chk ccall(
               (:PetscStrArrayallocpy, $petsc_library),
               PetscErrorCode,
               (Cchar, Cchar),
               list, t,
              )


	return nothing
end 

"""
	PetscStrArrayDestroy(petsclib::PetscLibType,list::Cchar) 
Frees array of strings created with `PetscStrArrayallocpy()`.

Not Collective; No Fortran Support

Output Parameter:
- `list` - array of strings

Level: intermediate

-seealso: `PetscStrArrayallocpy()`

# External Links
$(_doc_external("Sys/PetscStrArrayDestroy"))
"""
function PetscStrArrayDestroy(petsclib::PetscLibType, list::Cchar) end

@for_petsc function PetscStrArrayDestroy(petsclib::$UnionPetscLib, list::Cchar )

    @chk ccall(
               (:PetscStrArrayDestroy, $petsc_library),
               PetscErrorCode,
               (Cchar,),
               list,
              )


	return nothing
end 

"""
	PetscStrNArrayallocpy(petsclib::PetscLibType,n::PetscInt, list::Cchar, t::Cchar) 
Allocates space to hold a copy of an array of strings then copies the strings

Not Collective; No Fortran Support

Input Parameters:
- `n`    - the number of string entries
- `list` - pointer to array of strings

Output Parameter:
- `t` - the copied array string

Level: intermediate

-seealso: `PetscStrallocpy()`, `PetscStrArrayallocpy()`, `PetscStrNArrayDestroy()`

# External Links
$(_doc_external("Sys/PetscStrNArrayallocpy"))
"""
function PetscStrNArrayallocpy(petsclib::PetscLibType, n::PetscInt, list::Cchar, t::Cchar) end

@for_petsc function PetscStrNArrayallocpy(petsclib::$UnionPetscLib, n::$PetscInt, list::Cchar, t::Cchar )

    @chk ccall(
               (:PetscStrNArrayallocpy, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Cchar, Cchar),
               n, list, t,
              )


	return nothing
end 

"""
	n::PetscInt = PetscStrNArrayDestroy(petsclib::PetscLibType,list::Cchar) 
Frees array of strings created with `PetscStrNArrayallocpy()`.

Not Collective; No Fortran Support

Output Parameters:
- `n`    - number of string entries
- `list` - array of strings

Level: intermediate

-seealso: `PetscStrNArrayallocpy()`, `PetscStrArrayallocpy()`

# External Links
$(_doc_external("Sys/PetscStrNArrayDestroy"))
"""
function PetscStrNArrayDestroy(petsclib::PetscLibType, list::Cchar) end

@for_petsc function PetscStrNArrayDestroy(petsclib::$UnionPetscLib, list::Cchar )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscStrNArrayDestroy, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Cchar),
               n_, list,
              )

	n = n_[]

	return n
end 

"""
	t::PetscBool = PetscStrcasecmp(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}) 
Returns true if the two strings are the same
except possibly for case.

Not Collective; No Fortran Support

Input Parameters:
- `a` - pointer to first string
- `b` - pointer to second string

Output Parameter:
- `t` - if the two strings are the same

Level: intermediate

-seealso: `PetscStrcmp()`, `PetscStrncmp()`, `PetscStrgrt()`

# External Links
$(_doc_external("Sys/PetscStrcasecmp"))
"""
function PetscStrcasecmp(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}) end

@for_petsc function PetscStrcasecmp(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar} )
	t_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrcasecmp, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               a, b, t_,
              )

	t = t_[]

	return t
end 

"""
	cnt::PetscInt = PetscStrendswithwhich(petsclib::PetscLibType,a::Vector{Cchar}, bs::Cchar) 
Determines if a string ends with one of several possible strings

Not Collective; No Fortran Support

Input Parameters:
- `a`  - pointer to string
- `bs` - strings to end with (last entry must be `NULL`)

Output Parameter:
- `cnt` - the index of the string it ends with or the index of `NULL`

Level: intermediate

-seealso: `PetscStrbeginswithwhich()`, `PetscStrendswith()`, `PetscStrtoupper`, `PetscStrtolower()`, `PetscStrrchr()`, `PetscStrchr()`,
`PetscStrncmp()`, `PetscStrlen()`, `PetscStrcmp()`

# External Links
$(_doc_external("Sys/PetscStrendswithwhich"))
"""
function PetscStrendswithwhich(petsclib::PetscLibType, a::Vector{Cchar}, bs::Cchar) end

@for_petsc function PetscStrendswithwhich(petsclib::$UnionPetscLib, a::Vector{Cchar}, bs::Cchar )
	cnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscStrendswithwhich, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{$PetscInt}),
               a, bs, cnt_,
              )

	cnt = cnt_[]

	return cnt
end 

"""
	found::PetscBool = PetscStrInList(petsclib::PetscLibType,str::Vector{Cchar}, list::Vector{Cchar}, sep::Cchar) 
search for a string in character

Not Collective; No Fortran Support

Input Parameters:
- `str`  - the string to look for
- `list` - the list to search in
- `sep`  - the separator character

Output Parameter:
- `found` - whether `str` is in `list`

Level: intermediate

-seealso: `PetscTokenCreate()`, `PetscTokenFind()`, `PetscStrcmp()`

# External Links
$(_doc_external("Sys/PetscStrInList"))
"""
function PetscStrInList(petsclib::PetscLibType, str::Vector{Cchar}, list::Vector{Cchar}, sep::Cchar) end

@for_petsc function PetscStrInList(petsclib::$UnionPetscLib, str::Vector{Cchar}, list::Vector{Cchar}, sep::Cchar )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrInList, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Cchar, Ptr{PetscBool}),
               str, list, sep, found_,
              )

	found = found_[]

	return found
end 

"""
	PetscGetPetscDir(petsclib::PetscLibType,dir::Vector{Cchar}) 
Gets the directory PETSc is installed in

Not Collective; No Fortran Support

Output Parameter:
- `dir` - the directory

Level: developer

-seealso: `PetscGetArchType()`

# External Links
$(_doc_external("Sys/PetscGetPetscDir"))
"""
function PetscGetPetscDir(petsclib::PetscLibType, dir::Vector{Cchar}) end

@for_petsc function PetscGetPetscDir(petsclib::$UnionPetscLib, dir::Vector{Cchar} )
	dir_ = Ref(pointer(dir))

    @chk ccall(
               (:PetscGetPetscDir, $petsc_library),
               PetscErrorCode,
               (Ptr{Ptr{Cchar}},),
               dir_,
              )


	return nothing
end 

"""
	PetscStrreplace(petsclib::PetscLibType,comm::MPI_Comm, aa::Vector{Cchar}, b::Vector{Cchar}, len::Csize_t) 
Replaces substrings in string with other substrings

Not Collective; No Fortran Support

Input Parameters:
- `comm` - `MPI_Comm` of processors that are processing the string
- `aa`   - the string to look in
- `b`    - the resulting copy of a with replaced strings (`b` can be the same as `a`)
- `len`  - the length of `b`

Level: developer

-seealso: `PetscStrcmp()`

# External Links
$(_doc_external("Sys/PetscStrreplace"))
"""
function PetscStrreplace(petsclib::PetscLibType, comm::MPI_Comm, aa::Vector{Cchar}, b::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscStrreplace(petsclib::$UnionPetscLib, comm::MPI_Comm, aa::Vector{Cchar}, b::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscStrreplace, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               comm, aa, b, len,
              )


	return nothing
end 

"""
	value::PetscInt,found::PetscBool = PetscEListFind(petsclib::PetscLibType,n::PetscInt, list::Cchar, str::Vector{Cchar}) 
searches list of strings for given string, using case insensitive matching

Not Collective; No Fortran Support

Input Parameters:
- `n`    - number of strings in
- `list` - list of strings to search
- `str`  - string to look for, empty string "" accepts default (first entry in list)

Output Parameters:
- `value` - index of matching string (if found)
- `found` - boolean indicating whether string was found (can be `NULL`)

Level: developer

-seealso: `PetscEnumFind()`

# External Links
$(_doc_external("Sys/PetscEListFind"))
"""
function PetscEListFind(petsclib::PetscLibType, n::PetscInt, list::Cchar, str::Vector{Cchar}) end

@for_petsc function PetscEListFind(petsclib::$UnionPetscLib, n::$PetscInt, list::Cchar, str::Vector{Cchar} )
	value_ = Ref{$PetscInt}()
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscEListFind, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Cchar, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, list, str, value_, found_,
              )

	value = value_[]
	found = found_[]

	return value,found
end 

"""
	found::PetscBool = PetscEnumFind(petsclib::PetscLibType,enumlist::Cchar, str::Vector{Cchar}, value::PetscEnum) 
searches enum list of strings for given string, using case insensitive matching

Not Collective; No Fortran Support

Input Parameters:
- `enumlist` - list of strings to search, followed by enum name, then enum prefix, then `NULL`
- `str`      - string to look for

Output Parameters:
- `value` - index of matching string (if found)
- `found` - boolean indicating whether string was found (can be `NULL`)

Level: advanced

-seealso: `PetscEListFind()`

# External Links
$(_doc_external("Sys/PetscEnumFind"))
"""
function PetscEnumFind(petsclib::PetscLibType, enumlist::Cchar, str::Vector{Cchar}, value::PetscEnum) end

@for_petsc function PetscEnumFind(petsclib::$UnionPetscLib, enumlist::Cchar, str::Vector{Cchar}, value::PetscEnum )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscEnumFind, $petsc_library),
               PetscErrorCode,
               (Cchar, Ptr{Cchar}, Ptr{PetscEnum}, Ptr{PetscBool}),
               enumlist, str, value, found_,
              )

	found = found_[]

	return found
end 

"""
	PetscStrcat(petsclib::PetscLibType,s::Vector{Cchar}, t::Vector{Cchar}) 
Concatenates a string onto a given string

Not Collective, No Fortran Support

Input Parameters:
- `s` - string to be added to
- `t` - pointer to string to be added to end

Level: deprecated (since 3.18.5)

-seealso: `PetscStrlcat()`

# External Links
$(_doc_external("Sys/PetscStrcat"))
"""
function PetscStrcat(petsclib::PetscLibType, s::Vector{Cchar}, t::Vector{Cchar}) end

@for_petsc function PetscStrcat(petsclib::$UnionPetscLib, s::Vector{Cchar}, t::Vector{Cchar} )

    @chk ccall(
               (:PetscStrcat, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               s, t,
              )


	return nothing
end 

"""
	PetscStrcpy(petsclib::PetscLibType,s::Vector{Cchar}, t::Vector{Cchar}) 
Copies a string

Not Collective, No Fortran Support

Input Parameter:
- `t` - pointer to string

Output Parameter:
- `s` - the copied string

Level: deprecated (since 3.18.5)

-seealso: `PetscStrncpy()`

# External Links
$(_doc_external("Sys/PetscStrcpy"))
"""
function PetscStrcpy(petsclib::PetscLibType, s::Vector{Cchar}, t::Vector{Cchar}) end

@for_petsc function PetscStrcpy(petsclib::$UnionPetscLib, s::Vector{Cchar}, t::Vector{Cchar} )

    @chk ccall(
               (:PetscStrcpy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               s, t,
              )


	return nothing
end 

"""
	str::Vector{Cchar} = PetscGetArchType(petsclib::PetscLibType,slen::Csize_t) 
Returns the PETSC_ARCH that was used for this configuration of PETSc

Not Collective

Input Parameter:
- `slen` - length of string buffer

Output Parameter:
- `str` - string area to contain architecture name, should be at least 10 characters long. Name is truncated if string is not long enough.

Level: developer

-seealso: `PetscGetUserName()`, `PetscGetHostName()`

# External Links
$(_doc_external("Sys/PetscGetArchType"))
"""
function PetscGetArchType(petsclib::PetscLibType, slen::Csize_t) end

@for_petsc function PetscGetArchType(petsclib::$UnionPetscLib, slen::Csize_t )
	str = Vector{Cchar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscGetArchType, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               str, slen,
              )


	return str
end 

"""
	PetscSetDisplay(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscSetDisplay"))
"""
function PetscSetDisplay(petsclib::PetscLibType) end

@for_petsc function PetscSetDisplay(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSetDisplay, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscGetDisplay(petsclib::PetscLibType,display::Vector{Cchar}, n::Csize_t) 
Gets the X windows display variable for all processors.

Input Parameter:
- `n` - length of string display

Output Parameter:
- `display` - the display string

Options Database Keys:
- `-display <display>` - sets the display to use
- `-x_virtual`         - forces use of a X virtual display Xvfb that will not display anything but -draw_save will still work. Xvfb is automatically
started up in PetscSetDisplay() with this option

Level: advanced

-seealso: `PETSC_DRAW_X`, `PetscDrawOpenX()`

# External Links
$(_doc_external("Sys/PetscGetDisplay"))
"""
function PetscGetDisplay(petsclib::PetscLibType, display::Vector{Cchar}, n::Csize_t) end

@for_petsc function PetscGetDisplay(petsclib::$UnionPetscLib, display::Vector{Cchar}, n::Csize_t )

    @chk ccall(
               (:PetscGetDisplay, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               display, n,
              )


	return nothing
end 

"""
	PetscGetHostName(petsclib::PetscLibType,name::Vector{Cchar}, nlen::Csize_t) 
Returns the name of the host. This attempts to
return the entire Internet name. It may not return the same name
as `MPI_Get_processor_name()`.

Not Collective

Input Parameter:
- `nlen` - length of name

Output Parameter:
- `name` - contains host name.  Must be long enough to hold the name
This is the fully qualified name, including the domain.

Level: developer

-seealso: `PetscGetUserName()`, `PetscGetArchType()`

# External Links
$(_doc_external("Sys/PetscGetHostName"))
"""
function PetscGetHostName(petsclib::PetscLibType, name::Vector{Cchar}, nlen::Csize_t) end

@for_petsc function PetscGetHostName(petsclib::$UnionPetscLib, name::Vector{Cchar}, nlen::Csize_t )

    @chk ccall(
               (:PetscGetHostName, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               name, nlen,
              )


	return nothing
end 

"""
	sorted::PetscBool = PetscSortedReal(petsclib::PetscLibType,n::PetscCount, X::Vector{PetscReal}) 
Determines whether the array of `PetscReal` is sorted.

Not Collective

Input Parameters:
- `n` - number of values
- `X` - array of integers

Output Parameter:
- `sorted` - flag whether the array is sorted

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortedInt()`, `PetscSortedMPIInt()`

# External Links
$(_doc_external("Sys/PetscSortedReal"))
"""
function PetscSortedReal(petsclib::PetscLibType, n::PetscCount, X::Vector{PetscReal}) end

@for_petsc function PetscSortedReal(petsclib::$UnionPetscLib, n::PetscCount, X::Vector{$PetscReal} )
	sorted_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSortedReal, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscReal}, Ptr{PetscBool}),
               n, X, sorted_,
              )

	sorted = sorted_[]

	return sorted
end 

"""
	PetscSortReal(petsclib::PetscLibType,n::PetscCount, v::Vector{PetscReal}) 
Sorts an array of `PetscReal` in place in increasing order.

Not Collective

Input Parameters:
- `n` - number of values
- `v` - array of doubles

Level: intermediate

-seealso: `PetscRealSortSemiOrdered()`, `PetscSortInt()`, `PetscSortRealWithPermutation()`, `PetscSortRealWithArrayInt()`

# External Links
$(_doc_external("Sys/PetscSortReal"))
"""
function PetscSortReal(petsclib::PetscLibType, n::PetscCount, v::Vector{PetscReal}) end

@for_petsc function PetscSortReal(petsclib::$UnionPetscLib, n::PetscCount, v::Vector{$PetscReal} )

    @chk ccall(
               (:PetscSortReal, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscReal}),
               n, v,
              )


	return nothing
end 

"""
	PetscSortRealWithArrayInt(petsclib::PetscLibType,n::PetscCount, r::Vector{PetscReal}, Ii::Vector{PetscInt}) 
Sorts an array of `PetscReal` in place in increasing order;
changes a second `PetscInt` array to match the sorted first array.

Not Collective

Input Parameters:
- `n`  - number of values
- `Ii` - array of integers
- `r`  - second array of integers

Level: intermediate

-seealso: `PetscSortReal()`

# External Links
$(_doc_external("Sys/PetscSortRealWithArrayInt"))
"""
function PetscSortRealWithArrayInt(petsclib::PetscLibType, n::PetscCount, r::Vector{PetscReal}, Ii::Vector{PetscInt}) end

@for_petsc function PetscSortRealWithArrayInt(petsclib::$UnionPetscLib, n::PetscCount, r::Vector{$PetscReal}, Ii::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortRealWithArrayInt, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{$PetscReal}, Ptr{$PetscInt}),
               n, r, Ii,
              )


	return nothing
end 

"""
	loc::PetscInt = PetscFindReal(petsclib::PetscLibType,key::PetscReal, n::PetscCount, t::Vector{PetscReal}, eps::PetscReal) 
Finds a `PetscReal` in a sorted array of `PetscReal`s

Not Collective

Input Parameters:
- `key` - the value to locate
- `n`   - number of values in the array
- `t`   - array of values
- `eps` - tolerance used to compare

Output Parameter:
- `loc` - the location if found, otherwise -(slot+1) where slot is the place the value would go

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortRealWithArrayInt()`

# External Links
$(_doc_external("Sys/PetscFindReal"))
"""
function PetscFindReal(petsclib::PetscLibType, key::PetscReal, n::PetscCount, t::Vector{PetscReal}, eps::PetscReal) end

@for_petsc function PetscFindReal(petsclib::$UnionPetscLib, key::$PetscReal, n::PetscCount, t::Vector{$PetscReal}, eps::$PetscReal )
	loc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscFindReal, $petsc_library),
               PetscErrorCode,
               ($PetscReal, PetscCount, Ptr{$PetscReal}, $PetscReal, Ptr{$PetscInt}),
               key, n, t, eps, loc_,
              )

	loc = loc_[]

	return loc
end 

"""
	PetscSortRemoveDupsReal(petsclib::PetscLibType,n::PetscInt, v::Vector{PetscReal}) 
Sorts an array of `PetscReal` in place in increasing order and removes all duplicate entries

Not Collective

Input Parameters:
- `n` - initial number of values
- `v` - array of values

-seealso: `PetscSortReal()`, `PetscSortRemoveDupsInt()`

# External Links
$(_doc_external("Sys/PetscSortRemoveDupsReal"))
"""
function PetscSortRemoveDupsReal(petsclib::PetscLibType, n::PetscInt, v::Vector{PetscReal}) end

@for_petsc function PetscSortRemoveDupsReal(petsclib::$UnionPetscLib, n::$PetscInt, v::Vector{$PetscReal} )

    @chk ccall(
               (:PetscSortRemoveDupsReal, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscInt}, Ptr{$PetscReal}),
               n, v,
              )


	return nothing
end 

"""
	PetscSortSplit(petsclib::PetscLibType,ncut::PetscInt, n::PetscInt, a::Vector{PetscScalar}, idx::Vector{PetscInt}) 
Quick

Not Collective

Input Parameters:
- `ncut` - splitting index
- `n`    - number of values to sort

Input/Output Parameters:
- `a`   - array of values, on output the values are permuted such that its elements satisfy:
abs(a[i]) >= abs(a[ncut-1]) for i < ncut and
abs(a[i]) <= abs(a[ncut-1]) for i >= ncut
- `idx` - index for array a, on output permuted accordingly

Level: intermediate

-seealso: `PetscSortInt()`, `PetscSortRealWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortSplit"))
"""
function PetscSortSplit(petsclib::PetscLibType, ncut::PetscInt, n::PetscInt, a::Vector{PetscScalar}, idx::Vector{PetscInt}) end

@for_petsc function PetscSortSplit(petsclib::$UnionPetscLib, ncut::$PetscInt, n::$PetscInt, a::Vector{$PetscScalar}, idx::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortSplit, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscInt}),
               ncut, n, a, idx,
              )


	return nothing
end 

"""
	PetscSortSplitReal(petsclib::PetscLibType,ncut::PetscInt, n::PetscInt, a::Vector{PetscReal}, idx::Vector{PetscInt}) 
Quick

Not Collective

Input Parameters:
- `ncut` - splitting index
- `n`    - number of values to sort

Input/Output Parameters:
- `a`   - array of values, on output the values are permuted such that its elements satisfy:
abs(a[i]) >= abs(a[ncut-1]) for i < ncut and
abs(a[i]) <= abs(a[ncut-1]) for i >= ncut
- `idx` - index for array a, on output permuted accordingly

Level: intermediate

-seealso: `PetscSortInt()`, `PetscSortRealWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortSplitReal"))
"""
function PetscSortSplitReal(petsclib::PetscLibType, ncut::PetscInt, n::PetscInt, a::Vector{PetscReal}, idx::Vector{PetscInt}) end

@for_petsc function PetscSortSplitReal(petsclib::$UnionPetscLib, ncut::$PetscInt, n::$PetscInt, a::Vector{$PetscReal}, idx::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortSplitReal, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
               ncut, n, a, idx,
              )


	return nothing
end 

"""
	PetscGatherNumberOfMessages(petsclib::PetscLibType,comm::MPI_Comm, iflags::Vector{PetscMPIInt}, ilengths::Vector{PetscMPIInt}, nrecvs::PetscMPIInt) 
Computes the number of messages an MPI rank expects to receive during a neighbor communication

Collective, No Fortran Support

Input Parameters:
- `comm`     - Communicator
- `iflags`   - an array of integers of length sizeof(comm). A '1' in `ilengths`[i] represent a
message from current node to ith node. Optionally `NULL`
- `ilengths` - Non zero ilengths[i] represent a message to i of length `ilengths`[i].
Optionally `NULL`.

Output Parameter:
- `nrecvs` - number of messages received

Level: developer

-seealso: `PetscGatherMessageLengths()`, `PetscGatherMessageLengths2()`, `PetscCommBuildTwoSided()`

# External Links
$(_doc_external("Sys/PetscGatherNumberOfMessages"))
"""
function PetscGatherNumberOfMessages(petsclib::PetscLibType, comm::MPI_Comm, iflags::Vector{PetscMPIInt}, ilengths::Vector{PetscMPIInt}, nrecvs::PetscMPIInt) end

@for_petsc function PetscGatherNumberOfMessages(petsclib::$UnionPetscLib, comm::MPI_Comm, iflags::Vector{PetscMPIInt}, ilengths::Vector{PetscMPIInt}, nrecvs::PetscMPIInt )

    @chk ccall(
               (:PetscGatherNumberOfMessages, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}),
               comm, iflags, ilengths, nrecvs,
              )


	return nothing
end 

"""
	PetscGatherMessageLengths(petsclib::PetscLibType,comm::MPI_Comm, nsends::PetscMPIInt, nrecvs::PetscMPIInt, ilengths::Vector{PetscMPIInt}, onodes::PetscMPIInt, olengths::PetscMPIInt) 
Computes information about messages that an MPI rank will receive,
including (from-id,length) pairs for each message.

Collective, No Fortran Support

Input Parameters:
- `comm`     - Communicator
- `nsends`   - number of messages that are to be sent.
- `nrecvs`   - number of messages being received
- `ilengths` - an array of integers of length sizeof(comm)
a non zero `ilengths`[i] represent a message to i of length `ilengths`[i]

Output Parameters:
- `onodes`   - list of ranks from which messages are expected
- `olengths` - corresponding message lengths

Level: developer

-seealso: `PetscGatherNumberOfMessages()`, `PetscGatherMessageLengths2()`, `PetscCommBuildTwoSided()`

# External Links
$(_doc_external("Sys/PetscGatherMessageLengths"))
"""
function PetscGatherMessageLengths(petsclib::PetscLibType, comm::MPI_Comm, nsends::PetscMPIInt, nrecvs::PetscMPIInt, ilengths::Vector{PetscMPIInt}, onodes::PetscMPIInt, olengths::PetscMPIInt) end

@for_petsc function PetscGatherMessageLengths(petsclib::$UnionPetscLib, comm::MPI_Comm, nsends::PetscMPIInt, nrecvs::PetscMPIInt, ilengths::Vector{PetscMPIInt}, onodes::PetscMPIInt, olengths::PetscMPIInt )

    @chk ccall(
               (:PetscGatherMessageLengths, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscMPIInt, PetscMPIInt, Ptr{PetscMPIInt}, PetscMPIInt, PetscMPIInt),
               comm, nsends, nrecvs, ilengths, onodes, olengths,
              )


	return nothing
end 

"""
	PetscGatherMessageLengths2(petsclib::PetscLibType,comm::MPI_Comm, nsends::PetscMPIInt, nrecvs::PetscMPIInt, ilengths1::Vector{PetscMPIInt}, ilengths2::Vector{PetscMPIInt}, onodes::PetscMPIInt, olengths1::PetscMPIInt, olengths2::PetscMPIInt) 
Computes info about messages that a MPI rank will receive,
including (from-id,length) pairs for each message. Same functionality as `PetscGatherMessageLengths()`
except it takes TWO ilenths and output TWO olengths.

Collective, No Fortran Support

Input Parameters:
- `comm`      - Communicator
- `nsends`    - number of messages that are to be sent.
- `nrecvs`    - number of messages being received
- `ilengths1` - first array of integers of length sizeof(comm)
- `ilengths2` - second array of integers of length sizeof(comm)

Output Parameters:
- `onodes`    - list of ranks from which messages are expected
- `olengths1` - first corresponding message lengths
- `olengths2` - second  message lengths

Level: developer

-seealso: `PetscGatherMessageLengths()`, `PetscGatherNumberOfMessages()`, `PetscCommBuildTwoSided()`

# External Links
$(_doc_external("Sys/PetscGatherMessageLengths2"))
"""
function PetscGatherMessageLengths2(petsclib::PetscLibType, comm::MPI_Comm, nsends::PetscMPIInt, nrecvs::PetscMPIInt, ilengths1::Vector{PetscMPIInt}, ilengths2::Vector{PetscMPIInt}, onodes::PetscMPIInt, olengths1::PetscMPIInt, olengths2::PetscMPIInt) end

@for_petsc function PetscGatherMessageLengths2(petsclib::$UnionPetscLib, comm::MPI_Comm, nsends::PetscMPIInt, nrecvs::PetscMPIInt, ilengths1::Vector{PetscMPIInt}, ilengths2::Vector{PetscMPIInt}, onodes::PetscMPIInt, olengths1::PetscMPIInt, olengths2::PetscMPIInt )

    @chk ccall(
               (:PetscGatherMessageLengths2, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscMPIInt, PetscMPIInt, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}, PetscMPIInt, PetscMPIInt, PetscMPIInt),
               comm, nsends, nrecvs, ilengths1, ilengths2, onodes, olengths1, olengths2,
              )


	return nothing
end 

"""
	PetscPostIrecvInt(petsclib::PetscLibType,comm::MPI_Comm, tag::PetscMPIInt, nrecvs::PetscMPIInt, onodes::Vector{PetscMPIInt}, olengths::Vector{PetscMPIInt}, rbuf::PetscInt, r_waits::MPI_Request) 

# External Links
$(_doc_external("Sys/PetscPostIrecvInt"))
"""
function PetscPostIrecvInt(petsclib::PetscLibType, comm::MPI_Comm, tag::PetscMPIInt, nrecvs::PetscMPIInt, onodes::Vector{PetscMPIInt}, olengths::Vector{PetscMPIInt}, rbuf::PetscInt, r_waits::MPI_Request) end

@for_petsc function PetscPostIrecvInt(petsclib::$UnionPetscLib, comm::MPI_Comm, tag::PetscMPIInt, nrecvs::PetscMPIInt, onodes::Vector{PetscMPIInt}, olengths::Vector{PetscMPIInt}, rbuf::$PetscInt, r_waits::MPI_Request )

    @chk ccall(
               (:PetscPostIrecvInt, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscMPIInt, PetscMPIInt, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}, $PetscInt, MPI_Request),
               comm, tag, nrecvs, onodes, olengths, rbuf, r_waits,
              )


	return nothing
end 

"""
	PetscPostIrecvScalar(petsclib::PetscLibType,comm::MPI_Comm, tag::PetscMPIInt, nrecvs::PetscMPIInt, onodes::Vector{PetscMPIInt}, olengths::Vector{PetscMPIInt}, rbuf::PetscScalar, r_waits::MPI_Request) 

# External Links
$(_doc_external("Sys/PetscPostIrecvScalar"))
"""
function PetscPostIrecvScalar(petsclib::PetscLibType, comm::MPI_Comm, tag::PetscMPIInt, nrecvs::PetscMPIInt, onodes::Vector{PetscMPIInt}, olengths::Vector{PetscMPIInt}, rbuf::PetscScalar, r_waits::MPI_Request) end

@for_petsc function PetscPostIrecvScalar(petsclib::$UnionPetscLib, comm::MPI_Comm, tag::PetscMPIInt, nrecvs::PetscMPIInt, onodes::Vector{PetscMPIInt}, olengths::Vector{PetscMPIInt}, rbuf::$PetscScalar, r_waits::MPI_Request )

    @chk ccall(
               (:PetscPostIrecvScalar, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscMPIInt, PetscMPIInt, Ptr{PetscMPIInt}, Ptr{PetscMPIInt}, $PetscScalar, MPI_Request),
               comm, tag, nrecvs, onodes, olengths, rbuf, r_waits,
              )


	return nothing
end 

"""
	PetscShmgetAddressesFinalize(petsclib::PetscLibType) 
frees any shared memory that was allocated by `PetscShmgetAllocateArray()` but
not deallocated with `PetscShmgetDeallocateArray()`

Level: developer

-seealso: `PetscShmgetAllocateArray()`, `PetscShmgetDeallocateArray()`, `PetscShmgetUnmapAddresses()`

# External Links
$(_doc_external("Sys/PetscShmgetAddressesFinalize"))
"""
function PetscShmgetAddressesFinalize(petsclib::PetscLibType) end

@for_petsc function PetscShmgetAddressesFinalize(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscShmgetAddressesFinalize, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscShmgetMapAddresses(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, baseaddres::Cvoid, addres::Cvoid) 
given shared address on the first MPI process determines the
addresses on the other MPI processes that map to the same physical memory

Input Parameters:
- `comm`       - the `MPI_Comm` to scatter the address
- `n`          - the number of addresses, each obtained on MPI process zero by `PetscShmgetAllocateArray()`
- `baseaddres` - the addresses on the first MPI process, ignored on all but first process

Output Parameter:
- `addres` - the addresses on each MPI process, the array of void * must already be allocated

Level: developer

-seealso: `PetscShmgetDeallocateArray()`, `PetscShmgetAllocateArray()`, `PetscShmgetUnmapAddresses()`

# External Links
$(_doc_external("Sys/PetscShmgetMapAddresses"))
"""
function PetscShmgetMapAddresses(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, baseaddres::Cvoid, addres::Cvoid) end

@for_petsc function PetscShmgetMapAddresses(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, baseaddres::Cvoid, addres::Cvoid )

    @chk ccall(
               (:PetscShmgetMapAddresses, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Cvoid, Cvoid),
               comm, n, baseaddres, addres,
              )


	return nothing
end 

"""
	PetscShmgetUnmapAddresses(petsclib::PetscLibType,n::PetscInt, addres::Cvoid) 
given shared addresses on a MPI process unlink it

Input Parameters:
- `n`      - the number of addresses, each obtained on MPI process zero by `PetscShmgetAllocateArray()`
- `addres` - the addresses

Level: developer

-seealso: `PetscShmgetDeallocateArray()`, `PetscShmgetAllocateArray()`, `PetscShmgetMapAddresses()`

# External Links
$(_doc_external("Sys/PetscShmgetUnmapAddresses"))
"""
function PetscShmgetUnmapAddresses(petsclib::PetscLibType, n::PetscInt, addres::Cvoid) end

@for_petsc function PetscShmgetUnmapAddresses(petsclib::$UnionPetscLib, n::$PetscInt, addres::Cvoid )

    @chk ccall(
               (:PetscShmgetUnmapAddresses, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Cvoid),
               n, addres,
              )


	return nothing
end 

"""
	PetscShmgetAllocateArray(petsclib::PetscLibType,sz::Csize_t, asz::Csize_t, addr::Vector{Cvoid}) 
allocates shared memory accessible by all MPI processes in the server

Not Collective, only called on the first MPI process

Input Parameters:
- `sz`  - the number of elements in the array
- `asz` - the size of an entry in the array, for example `sizeof(PetscScalar)`

Output Parameters:
- `addr` - the address of the array

Level: developer

-seealso: [](sec_pcmpi), `PCMPIServerBegin()`, `PCMPI`, `KSPCheckPCMPI()`, `PetscShmgetDeallocateArray()`

# External Links
$(_doc_external("Sys/PetscShmgetAllocateArray"))
"""
function PetscShmgetAllocateArray(petsclib::PetscLibType, sz::Csize_t, asz::Csize_t, addr::Vector{Cvoid}) end

@for_petsc function PetscShmgetAllocateArray(petsclib::$UnionPetscLib, sz::Csize_t, asz::Csize_t, addr::Vector{Cvoid} )
	addr_ = Ref(pointer(addr))

    @chk ccall(
               (:PetscShmgetAllocateArray, $petsc_library),
               PetscErrorCode,
               (Csize_t, Csize_t, Ptr{Ptr{Cvoid}}),
               sz, asz, addr_,
              )


	return nothing
end 

"""
	PetscShmgetDeallocateArray(petsclib::PetscLibType,addr::Vector{Cvoid}) 
deallocates shared memory accessible by all MPI processes in the server

Not Collective, only called on the first MPI process

Input Parameter:
- `addr` - the address of array

Level: developer

-seealso: [](sec_pcmpi), `PCMPIServerBegin()`, `PCMPI`, `KSPCheckPCMPI()`, `PetscShmgetAllocateArray()`

# External Links
$(_doc_external("Sys/PetscShmgetDeallocateArray"))
"""
function PetscShmgetDeallocateArray(petsclib::PetscLibType, addr::Vector{Cvoid}) end

@for_petsc function PetscShmgetDeallocateArray(petsclib::$UnionPetscLib, addr::Vector{Cvoid} )
	addr_ = Ref(pointer(addr))

    @chk ccall(
               (:PetscShmgetDeallocateArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Ptr{Cvoid}},),
               addr_,
              )


	return nothing
end 

"""
	PetscMPIDump(petsclib::PetscLibType,fd::Libc.FILE) 
Dumps a listing of incomplete MPI operations, such as sends that
have never been received, etc.

Collective on `PETSC_COMM_WORLD`

Input Parameter:
- `fd` - file pointer.  If fp is `NULL`, `stdout` is assumed.

Options Database Key:
- `-mpidump` - Dumps MPI incompleteness during call to PetscFinalize()

Level: developer

-seealso: `PetscMallocDump()`

# External Links
$(_doc_external("Sys/PetscMPIDump"))
"""
function PetscMPIDump(petsclib::PetscLibType, fd::Libc.FILE) end

@for_petsc function PetscMPIDump(petsclib::$UnionPetscLib, fd::Libc.FILE )

    @chk ccall(
               (:PetscMPIDump, $petsc_library),
               PetscErrorCode,
               (Ptr{Libc.FILE},),
               fd,
              )


	return nothing
end 

"""
	PetscSortIntWithPermutation(petsclib::PetscLibType,n::PetscInt, i::Vector{PetscInt}, idx::Vector{PetscInt}) 
Computes the permutation of `PetscInt` that gives
a sorted sequence.

Not Collective

Input Parameters:
- `n`   - number of values to sort
- `i`   - values to sort
- `idx` - permutation array. Must be initialized to 0:`n`-1 on input.

Level: intermediate

-seealso: `PetscSortInt()`, `PetscSortRealWithPermutation()`, `PetscSortIntWithArray()`

# External Links
$(_doc_external("Sys/PetscSortIntWithPermutation"))
"""
function PetscSortIntWithPermutation(petsclib::PetscLibType, n::PetscInt, i::Vector{PetscInt}, idx::Vector{PetscInt}) end

@for_petsc function PetscSortIntWithPermutation(petsclib::$UnionPetscLib, n::$PetscInt, i::Vector{$PetscInt}, idx::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortIntWithPermutation, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, i, idx,
              )


	return nothing
end 

"""
	PetscSortRealWithPermutation(petsclib::PetscLibType,n::PetscInt, i::Vector{PetscReal}, idx::Vector{PetscInt}) 
Computes the permutation of `PetscReal` that gives
a sorted sequence.

Not Collective

Input Parameters:
- `n`   - number of values to sort
- `i`   - values to sort
- `idx` - permutation array. Must be initialized to 0:`n`-1 on input.

Level: intermediate

-seealso: `PetscSortReal()`, `PetscSortIntWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortRealWithPermutation"))
"""
function PetscSortRealWithPermutation(petsclib::PetscLibType, n::PetscInt, i::Vector{PetscReal}, idx::Vector{PetscInt}) end

@for_petsc function PetscSortRealWithPermutation(petsclib::$UnionPetscLib, n::$PetscInt, i::Vector{$PetscReal}, idx::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSortRealWithPermutation, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscInt}),
               n, i, idx,
              )


	return nothing
end 

"""
	PetscSortStrWithPermutation(petsclib::PetscLibType,n::PetscInt, i::Vector{Cchar}, idx::Vector{PetscInt}) 
Computes the permutation of strings that gives
a sorted sequence.

Not Collective, No Fortran Support

Input Parameters:
- `n`   - number of values to sort
- `i`   - values to sort
- `idx` - permutation array. Must be initialized to 0:`n`-1 on input.

Level: intermediate

-seealso: `PetscSortInt()`, `PetscSortRealWithPermutation()`

# External Links
$(_doc_external("Sys/PetscSortStrWithPermutation"))
"""
function PetscSortStrWithPermutation(petsclib::PetscLibType, n::PetscInt, i::Vector{Cchar}, idx::Vector{PetscInt}) end

@for_petsc function PetscSortStrWithPermutation(petsclib::$UnionPetscLib, n::$PetscInt, i::Vector{Cchar}, idx::Vector{$PetscInt} )
	i_ = Ref(pointer(i))

    @chk ccall(
               (:PetscSortStrWithPermutation, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{Cchar}}, Ptr{$PetscInt}),
               n, i_, idx,
              )


	return nothing
end 

"""
	PetscSleep(petsclib::PetscLibType,s::PetscReal) 
Sleeps some number of seconds.

Not Collective

Input Parameter:
- `s` - number of seconds to sleep

Level: intermediate

-seealso: `PetscTime()`

# External Links
$(_doc_external("Sys/PetscSleep"))
"""
function PetscSleep(petsclib::PetscLibType, s::PetscReal) end

@for_petsc function PetscSleep(petsclib::$UnionPetscLib, s::$PetscReal )

    @chk ccall(
               (:PetscSleep, $petsc_library),
               PetscErrorCode,
               ($PetscReal,),
               s,
              )


	return nothing
end 

"""
	PetscBarrier(petsclib::PetscLibType,obj::PetscObject) 
Blocks until this routine is executed by all processors owning the object `obj`.

Input Parameter:
- `obj` - PETSc object  (`Mat`, `Vec`, `IS`, `SNES` etc...)

Level: intermediate

-seealso: `PetscObject`, `MPI_Comm`, `MPI_Barrier`

# External Links
$(_doc_external("Sys/PetscBarrier"))
"""
function PetscBarrier(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscBarrier(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscBarrier, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscSequentialPhaseBegin(petsclib::PetscLibType,comm::MPI_Comm, ng::Cint) 
Begins a sequential section of code.

Collective

Input Parameters:
- `comm` - Communicator to sequentialize over
- `ng`   - Number in processor group.  This many processes are allowed to execute
at the same time (usually 1)

Level: intermediate

-seealso: `PetscSequentialPhaseEnd()`, `PetscSynchronizedPrintf()`

# External Links
$(_doc_external("Sys/PetscSequentialPhaseBegin"))
"""
function PetscSequentialPhaseBegin(petsclib::PetscLibType, comm::MPI_Comm, ng::Cint) end

@for_petsc function PetscSequentialPhaseBegin(petsclib::$UnionPetscLib, comm::MPI_Comm, ng::Cint )

    @chk ccall(
               (:PetscSequentialPhaseBegin, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Cint),
               comm, ng,
              )


	return nothing
end 

"""
	PetscSequentialPhaseEnd(petsclib::PetscLibType,comm::MPI_Comm, ng::Cint) 
Ends a sequential section of code.

Collective

Input Parameters:
- `comm` - Communicator to sequentialize.
- `ng`   - Number in processor group.  This many processes are allowed to execute
at the same time (usually 1)

Level: intermediate

-seealso: `PetscSequentialPhaseBegin()`

# External Links
$(_doc_external("Sys/PetscSequentialPhaseEnd"))
"""
function PetscSequentialPhaseEnd(petsclib::PetscLibType, comm::MPI_Comm, ng::Cint) end

@for_petsc function PetscSequentialPhaseEnd(petsclib::$UnionPetscLib, comm::MPI_Comm, ng::Cint )

    @chk ccall(
               (:PetscSequentialPhaseEnd, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Cint),
               comm, ng,
              )


	return nothing
end 

"""
	minMaxValGlobal::Vector{PetscInt} = PetscGlobalMinMaxInt(petsclib::PetscLibType,comm::MPI_Comm, minMaxVal::Vector{PetscInt}) 
Get the global min/max from local min/max input

Collective

Input Parameters:
- `comm`      - The MPI communicator to reduce with
- `minMaxVal` - An array with the local min and max

Output Parameter:
- `minMaxValGlobal` - An array with the global min and max

Level: beginner

-seealso: `PetscSplitOwnership()`, `PetscGlobalMinMaxReal()`

# External Links
$(_doc_external("Sys/PetscGlobalMinMaxInt"))
"""
function PetscGlobalMinMaxInt(petsclib::PetscLibType, comm::MPI_Comm, minMaxVal::Vector{PetscInt}) end

@for_petsc function PetscGlobalMinMaxInt(petsclib::$UnionPetscLib, comm::MPI_Comm, minMaxVal::Vector{$PetscInt} )
	minMaxValGlobal = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscGlobalMinMaxInt, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{$PetscInt}, Ptr{$PetscInt}),
               comm, minMaxVal, minMaxValGlobal,
              )


	return minMaxValGlobal
end 

"""
	minMaxValGlobal::Vector{PetscReal} = PetscGlobalMinMaxReal(petsclib::PetscLibType,comm::MPI_Comm, minMaxVal::Vector{PetscReal}) 
Get the global min/max from local min/max input

Collective

Input Parameters:
- `comm`      - The MPI communicator to reduce with
- `minMaxVal` - An array with the local min and max

Output Parameter:
- `minMaxValGlobal` - An array with the global min and max

Level: beginner

-seealso: `PetscSplitOwnership()`, `PetscGlobalMinMaxInt()`

# External Links
$(_doc_external("Sys/PetscGlobalMinMaxReal"))
"""
function PetscGlobalMinMaxReal(petsclib::PetscLibType, comm::MPI_Comm, minMaxVal::Vector{PetscReal}) end

@for_petsc function PetscGlobalMinMaxReal(petsclib::$UnionPetscLib, comm::MPI_Comm, minMaxVal::Vector{$PetscReal} )
	minMaxValGlobal = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscGlobalMinMaxReal, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{$PetscReal}, Ptr{$PetscReal}),
               comm, minMaxVal, minMaxValGlobal,
              )


	return minMaxValGlobal
end 

"""
	e::PetscBool = PetscMemcmp(petsclib::PetscLibType,str1::Cvoid, str2::Cvoid, len::Csize_t) 
Compares two byte streams in memory.

Not Collective

Input Parameters:
- `str1` - Pointer to the first byte stream
- `str2` - Pointer to the second byte stream
- `len`  - The length of the byte stream
(both str1 and str2 are assumed to be of length len)

Output Parameter:
- `e` - `PETSC_TRUE` if equal else `PETSC_FALSE`.

Level: intermediate

-seealso: `PetscMemcpy()`, `PetscArrayzero()`, `PetscMemzero()`, `PetscArraycmp()`, `PetscArraycpy()`, `PetscStrallocpy()`,
`PetscArraymove()`

# External Links
$(_doc_external("Sys/PetscMemcmp"))
"""
function PetscMemcmp(petsclib::PetscLibType, str1::Cvoid, str2::Cvoid, len::Csize_t) end

@for_petsc function PetscMemcmp(petsclib::$UnionPetscLib, str1::Cvoid, str2::Cvoid, len::Csize_t )
	e_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscMemcmp, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t, Ptr{PetscBool}),
               str1, str2, len, e_,
              )

	e = e_[]

	return e
end 

"""
	PetscProcessPlacementView(petsclib::PetscLibType,viewer::PetscViewer) 
display the MPI rank placement by core

Input Parameter:
- `viewer` - `PETSCVIEWERASCII` to display the results on

Level: intermediate

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscProcessPlacementView"))
"""
function PetscProcessPlacementView(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscProcessPlacementView(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscProcessPlacementView, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscGetUserName(petsclib::PetscLibType,name::Vector{Cchar}, nlen::Csize_t) 
Returns the name of the user.

Not Collective

Input Parameter:
- `nlen` - length of name

Output Parameter:
- `name` - contains user name. Must be long enough to hold the name

Level: developer

-seealso: `PetscGetHostName()`

# External Links
$(_doc_external("Sys/PetscGetUserName"))
"""
function PetscGetUserName(petsclib::PetscLibType, name::Vector{Cchar}, nlen::Csize_t) end

@for_petsc function PetscGetUserName(petsclib::$UnionPetscLib, name::Vector{Cchar}, nlen::Csize_t )

    @chk ccall(
               (:PetscGetUserName, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               name, nlen,
              )


	return nothing
end 

"""
	PetscGetDate(petsclib::PetscLibType,date::Vector{Cchar}, len::Csize_t) 
Gets the current date.

Not Collective

Input Parameter:
- `len` - length of string to hold date

Output Parameter:
- `date` - the date

Level: beginner

-seealso: `PetscGetHostName()`

# External Links
$(_doc_external("Sys/PetscGetDate"))
"""
function PetscGetDate(petsclib::PetscLibType, date::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscGetDate(petsclib::$UnionPetscLib, date::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscGetDate, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Csize_t),
               date, len,
              )


	return nothing
end 

"""
	PetscGetCPUTime(petsclib::PetscLibType,t::PetscLogDouble) 
Returns the CPU time in seconds used by the process.

Not Collective

Output Parameter:
- `t`   - Time in seconds charged to the process.

Example:
-seealso: `PetscTime()`, `PetscLogView()`

# External Links
$(_doc_external("Sys/PetscGetCPUTime"))
"""
function PetscGetCPUTime(petsclib::PetscLibType, t::PetscLogDouble) end

@for_petsc function PetscGetCPUTime(petsclib::$UnionPetscLib, t::PetscLogDouble )

    @chk ccall(
               (:PetscGetCPUTime, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               t,
              )


	return nothing
end 

"""
	enabled::PetscBool = PetscInfoEnabled(petsclib::PetscLibType,classid::PetscClassId) 
Checks whether a given `PetscClassid` is allowed to print using `PetscInfo()`

Not Collective

Input Parameter:
- `classid` - `PetscClassid` retrieved from a `PetscObject` e.g. `VEC_CLASSID`

Output Parameter:
- `enabled` - `PetscBool` indicating whether this classid is allowed to print

Level: advanced

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoAllow()`, `PetscInfoGetInfo()`, `PetscObjectGetClassid()`

# External Links
$(_doc_external("Sys/PetscInfoEnabled"))
"""
function PetscInfoEnabled(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscInfoEnabled(petsclib::$UnionPetscLib, classid::PetscClassId )
	enabled_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscInfoEnabled, $petsc_library),
               PetscErrorCode,
               (PetscClassId, Ptr{PetscBool}),
               classid, enabled_,
              )

	enabled = enabled_[]

	return enabled
end 

"""
	PetscInfoAllow(petsclib::PetscLibType,flag::PetscBool) 
Enables/disables `PetscInfo()` messages

Not Collective

Input Parameter:
- `flag` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoEnabled()`, `PetscInfoGetInfo()`, `PetscInfoSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscInfoAllow"))
"""
function PetscInfoAllow(petsclib::PetscLibType, flag::PetscBool) end

@for_petsc function PetscInfoAllow(petsclib::$UnionPetscLib, flag::PetscBool )

    @chk ccall(
               (:PetscInfoAllow, $petsc_library),
               PetscErrorCode,
               (PetscBool,),
               flag,
              )


	return nothing
end 

"""
	PetscInfoSetFile(petsclib::PetscLibType,filename::Vector{Cchar}, mode::Vector{Cchar}) 
Sets the printing destination for all `PetscInfo()` calls

Not Collective

Input Parameters:
- `filename` - Name of the file where `PetscInfo()` will print to, use `NULL` to write to `PETSC_STDOUT`.
- `mode`     - Write mode passed to `PetscFOpen()`

Level: advanced

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoGetFile()`, `PetscInfoSetFromOptions()`, `PetscFOpen()`

# External Links
$(_doc_external("Sys/PetscInfoSetFile"))
"""
function PetscInfoSetFile(petsclib::PetscLibType, filename::Vector{Cchar}, mode::Vector{Cchar}) end

@for_petsc function PetscInfoSetFile(petsclib::$UnionPetscLib, filename::Vector{Cchar}, mode::Vector{Cchar} )

    @chk ccall(
               (:PetscInfoSetFile, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}),
               filename, mode,
              )


	return nothing
end 

"""
	PetscInfoGetFile(petsclib::PetscLibType,filename::Vector{Cchar}, InfoFile::Libc.FILE) 
Gets the `filename` and `FILE` pointer of the file where `PetscInfo()` prints to

Not Collective; No Fortran Support

Output Parameters:
- `filename` - The name of the output file
- `InfoFile` - The `FILE` pointer for the output file

Level: advanced

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoSetFile()`, `PetscInfoSetFromOptions()`, `PetscInfoDestroy()`

# External Links
$(_doc_external("Sys/PetscInfoGetFile"))
"""
function PetscInfoGetFile(petsclib::PetscLibType, filename::Vector{Cchar}, InfoFile::Libc.FILE) end

@for_petsc function PetscInfoGetFile(petsclib::$UnionPetscLib, filename::Vector{Cchar}, InfoFile::Libc.FILE )
	filename_ = Ref(pointer(filename))

    @chk ccall(
               (:PetscInfoGetFile, $petsc_library),
               PetscErrorCode,
               (Ptr{Ptr{Cchar}}, Libc.FILE),
               filename_, InfoFile,
              )


	return nothing
end 

"""
	PetscInfoSetClasses(petsclib::PetscLibType,exclude::PetscBool, n::PetscInt, classnames::Cchar) 
Sets the classes which `PetscInfo()` is filtered for/against

Not Collective; No Fortran Support

Input Parameters:
- `exclude`    - Whether or not to invert the filter, i.e. if exclude is true, `PetscInfo()` will print from every class that
is NOT one of the classes specified
- `n`          - Number of classes to filter for (size of `classnames`)
- `classnames` - String array containing the names of classes to filter for, e.g. "vec"

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoGetClass()`, `PetscInfoProcessClass()`, `PetscInfoSetFromOptions()`, `PetscStrToArray()`, `PetscObjectGetName()`

# External Links
$(_doc_external("Sys/PetscInfoSetClasses"))
"""
function PetscInfoSetClasses(petsclib::PetscLibType, exclude::PetscBool, n::PetscInt, classnames::Cchar) end

@for_petsc function PetscInfoSetClasses(petsclib::$UnionPetscLib, exclude::PetscBool, n::$PetscInt, classnames::Cchar )

    @chk ccall(
               (:PetscInfoSetClasses, $petsc_library),
               PetscErrorCode,
               (PetscBool, $PetscInt, Cchar),
               exclude, n, classnames,
              )


	return nothing
end 

"""
	found::PetscBool = PetscInfoGetClass(petsclib::PetscLibType,classname::Vector{Cchar}) 
Indicates whether the provided `classname` is marked as a filter in `PetscInfo()` as set by `PetscInfoSetClasses()`

Not Collective

Input Parameter:
- `classname` - Name of the class to search for

Output Parameter:
- `found` - `PetscBool` indicating whether the classname was found

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoSetClasses()`, `PetscInfoSetFromOptions()`, `PetscObjectGetName()`

# External Links
$(_doc_external("Sys/PetscInfoGetClass"))
"""
function PetscInfoGetClass(petsclib::PetscLibType, classname::Vector{Cchar}) end

@for_petsc function PetscInfoGetClass(petsclib::$UnionPetscLib, classname::Vector{Cchar} )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscInfoGetClass, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscBool}),
               classname, found_,
              )

	found = found_[]

	return found
end 

"""
	infoEnabled::PetscBool,classesSet::PetscBool,exclude::PetscBool,locked::PetscBool = PetscInfoGetInfo(petsclib::PetscLibType,commSelfFlag::PetscInfoCommFlag) 
Returns the current state of several flags for `PetscInfo()`

Not Collective

Output Parameters:
- `infoEnabled`  - `PETSC_TRUE` if `PetscInfoAllow`(`PETSC_TRUE`) has been called
- `classesSet`   - `PETSC_TRUE` if the list of classes to filter for has been set
- `exclude`      - `PETSC_TRUE` if the class filtering for `PetscInfo()` is inverted
- `locked`       - `PETSC_TRUE` if the list of classes to filter for has been locked
- `commSelfFlag` - Enum indicating whether `PetscInfo()` will print for communicators of size 1, any size != 1, or all
communicators

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoAllow()`, `PetscInfoSetFilterCommSelf`, `PetscInfoSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscInfoGetInfo"))
"""
function PetscInfoGetInfo(petsclib::PetscLibType, commSelfFlag::PetscInfoCommFlag) end

@for_petsc function PetscInfoGetInfo(petsclib::$UnionPetscLib, commSelfFlag::PetscInfoCommFlag )
	infoEnabled_ = Ref{PetscBool}()
	classesSet_ = Ref{PetscBool}()
	exclude_ = Ref{PetscBool}()
	locked_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscInfoGetInfo, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscInfoCommFlag}),
               infoEnabled_, classesSet_, exclude_, locked_, commSelfFlag,
              )

	infoEnabled = infoEnabled_[]
	classesSet = classesSet_[]
	exclude = exclude_[]
	locked = locked_[]

	return infoEnabled,classesSet,exclude,locked
end 

"""
	PetscInfoProcessClass(petsclib::PetscLibType,classname::Vector{Cchar}, numClassID::PetscInt, classIDs::Vector{PetscClassId}) 
Activates or deactivates a class based on the filtering status of `PetscInfo()`

Not Collective

Input Parameters:
- `classname`  - Name of the class to activate/deactivate `PetscInfo()` for
- `numClassID` - Number of entries in `classIDs`
- `classIDs`   - Array containing all of the `PetscClassId`s associated with `classname`

Options Database Key:
- `-info [filename][:[~]<list,of,classnames>[:[~]self]]` - specify which informative messages are printed, see `PetscInfo()`.

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoActivateClass()`, `PetscInfoDeactivateClass()`, `PetscInfoSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscInfoProcessClass"))
"""
function PetscInfoProcessClass(petsclib::PetscLibType, classname::Vector{Cchar}, numClassID::PetscInt, classIDs::Vector{PetscClassId}) end

@for_petsc function PetscInfoProcessClass(petsclib::$UnionPetscLib, classname::Vector{Cchar}, numClassID::$PetscInt, classIDs::Vector{PetscClassId} )

    @chk ccall(
               (:PetscInfoProcessClass, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, $PetscInt, Ptr{PetscClassId}),
               classname, numClassID, classIDs,
              )


	return nothing
end 

"""
	PetscInfoSetFilterCommSelf(petsclib::PetscLibType,commSelfFlag::PetscInfoCommFlag) 
Sets `PetscInfoCommFlag` enum to determine communicator filtering for `PetscInfo()`

Not Collective

Input Parameter:
- `commSelfFlag` - Enum value indicating method with which to filter `PetscInfo()` based on the size of the communicator of the object calling `PetscInfo()`

Options Database Key:
- `-info [filename][:[~]<list,of,classnames>[:[~]self]]` - specify which informative messages are printed, See `PetscInfo()`.

Level: advanced

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoGetInfo()`

# External Links
$(_doc_external("Sys/PetscInfoSetFilterCommSelf"))
"""
function PetscInfoSetFilterCommSelf(petsclib::PetscLibType, commSelfFlag::PetscInfoCommFlag) end

@for_petsc function PetscInfoSetFilterCommSelf(petsclib::$UnionPetscLib, commSelfFlag::PetscInfoCommFlag )

    @chk ccall(
               (:PetscInfoSetFilterCommSelf, $petsc_library),
               PetscErrorCode,
               (PetscInfoCommFlag,),
               commSelfFlag,
              )


	return nothing
end 

"""
	PetscInfoSetFromOptions(petsclib::PetscLibType,options::PetscOptions) 
Configure `PetscInfo()` using command line options, enabling or disabling various calls to `PetscInfo()`

Not Collective

Input Parameter:
- `options` - Options database, use `NULL` for default global database

Options Database Key:
- `-info [filename][:[~]<list,of,classnames>[:[~]self]]` - specify which informative messages are printed, See `PetscInfo()`.

Level: advanced

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoAllow()`, `PetscInfoSetFile()`, `PetscInfoSetClasses()`, `PetscInfoSetFilterCommSelf()`, `PetscInfoDestroy()`

# External Links
$(_doc_external("Sys/PetscInfoSetFromOptions"))
"""
function PetscInfoSetFromOptions(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscInfoSetFromOptions(petsclib::$UnionPetscLib, options::PetscOptions )

    @chk ccall(
               (:PetscInfoSetFromOptions, $petsc_library),
               PetscErrorCode,
               (COptions,),
               options,
              )

	return nothing
end 

"""
	PetscInfoDestroy(petsclib::PetscLibType) 
Destroys and resets internal `PetscInfo()` data structures.

Not Collective

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfo()`, `PetscInfoSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscInfoDestroy"))
"""
function PetscInfoDestroy(petsclib::PetscLibType) end

@for_petsc function PetscInfoDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscInfoDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscInfoDeactivateClass(petsclib::PetscLibType,classid::PetscClassId) 
Deactivates `PetscInfo()` messages for a PETSc object class.

Not Collective

Input Parameter:
- `classid` - The object class,  e.g., `MAT_CLASSID`, `SNES_CLASSID`, etc.

Options Database Key:
- `-info [filename][:[~]<list,of,classnames>[:[~]self]]` - specify which informative messages are printed, See `PetscInfo()`.

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfoActivateClass()`, `PetscInfo()`, `PetscInfoAllow()`, `PetscInfoSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscInfoDeactivateClass"))
"""
function PetscInfoDeactivateClass(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscInfoDeactivateClass(petsclib::$UnionPetscLib, classid::PetscClassId )

    @chk ccall(
               (:PetscInfoDeactivateClass, $petsc_library),
               PetscErrorCode,
               (PetscClassId,),
               classid,
              )


	return nothing
end 

"""
	PetscInfoActivateClass(petsclib::PetscLibType,classid::PetscClassId) 
Activates `PetscInfo()` messages for a PETSc object class.

Not Collective

Input Parameter:
- `classid` - The object class, e.g., `MAT_CLASSID`, `SNES_CLASSID`, etc.

Options Database Key:
- `-info [filename][:[~]<list,of,classnames>[:[~]self]]` - specify which informative messages are printed, See `PetscInfo()`.

Level: developer

-seealso: [](sec_PetscInfo), `PetscInfoDeactivateClass()`, `PetscInfo()`, `PetscInfoAllow()`, `PetscInfoSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscInfoActivateClass"))
"""
function PetscInfoActivateClass(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscInfoActivateClass(petsclib::$UnionPetscLib, classid::PetscClassId )

    @chk ccall(
               (:PetscInfoActivateClass, $petsc_library),
               PetscErrorCode,
               (PetscClassId,),
               classid,
              )


	return nothing
end 

"""
	PetscDLOpen(petsclib::PetscLibType,name::Vector{Cchar}, mode::PetscDLMode, handle::PetscDLHandle) 
opens a dynamic library

Not Collective, No Fortran Support

Input Parameters:
- `name` - name of library
- `mode` - options on how to open library

Output Parameter:
- `handle` - opaque pointer to be used with `PetscDLSym()`

Level: developer

-seealso: `PetscDLClose()`, `PetscDLSym()`, `PetscDLAddr()`, `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`,
`PetscDLLibraryRetrieve()`, `PetscDLLibraryOpen()`, `PetscDLLibraryClose()`, `PetscDLLibrarySym()`

# External Links
$(_doc_external("Sys/PetscDLOpen"))
"""
function PetscDLOpen(petsclib::PetscLibType, name::Vector{Cchar}, mode::PetscDLMode, handle::PetscDLHandle) end

@for_petsc function PetscDLOpen(petsclib::$UnionPetscLib, name::Vector{Cchar}, mode::PetscDLMode, handle::PetscDLHandle )

    @chk ccall(
               (:PetscDLOpen, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscDLMode, Ptr{PetscDLHandle}),
               name, mode, handle,
              )


	return nothing
end 

"""
	PetscDLClose(petsclib::PetscLibType,handle::PetscDLHandle) 
closes a dynamic library

Not Collective, No Fortran Support

Input Parameter:
- `handle` - the handle for the library obtained with `PetscDLOpen()`

Level: developer

-seealso: `PetscDLOpen()`, `PetscDLSym()`, `PetscDLAddr()`

# External Links
$(_doc_external("Sys/PetscDLClose"))
"""
function PetscDLClose(petsclib::PetscLibType, handle::PetscDLHandle) end

@for_petsc function PetscDLClose(petsclib::$UnionPetscLib, handle::PetscDLHandle )

    @chk ccall(
               (:PetscDLClose, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDLHandle},),
               handle,
              )


	return nothing
end 

"""
	PetscDLSym(petsclib::PetscLibType,handle::PetscDLHandle, symbol::Vector{Cchar}, value::Cvoid) 
finds a symbol in a dynamic library

Not Collective, No Fortran Support

Input Parameters:
- `handle` - obtained with `PetscDLOpen()` or `NULL`
- `symbol` - name of symbol

Output Parameter:
- `value` - pointer to the function, `NULL` if not found

Level: developer

-seealso: `PetscDLClose()`, `PetscDLOpen()`, `PetscDLAddr()`, `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`,
`PetscDLLibraryRetrieve()`, `PetscDLLibraryOpen()`, `PetscDLLibraryClose()`, `PetscDLLibrarySym()`

# External Links
$(_doc_external("Sys/PetscDLSym"))
"""
function PetscDLSym(petsclib::PetscLibType, handle::PetscDLHandle, symbol::Vector{Cchar}, value::Cvoid) end

@for_petsc function PetscDLSym(petsclib::$UnionPetscLib, handle::PetscDLHandle, symbol::Vector{Cchar}, value::Cvoid )

    @chk ccall(
               (:PetscDLSym, $petsc_library),
               PetscErrorCode,
               (PetscDLHandle, Ptr{Cchar}, Cvoid),
               handle, symbol, value,
              )


	return nothing
end 

"""
	PetscDLAddr(petsclib::PetscLibType,func::PetscVoidFn, name::Vector{Cchar}) 
find the name of a symbol in a dynamic library

Not Collective, No Fortran Support

Input Parameters:
- `func` - pointer to the function, `NULL` if not found

Output Parameter:
- `name` - name of symbol, or `NULL` if name lookup is not supported.

Level: developer

-seealso: `PetscDLClose()`, `PetscDLSym()`, `PetscDLOpen()`, `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`,
`PetscDLLibraryRetrieve()`, `PetscDLLibraryOpen()`, `PetscDLLibraryClose()`, `PetscDLLibrarySym()`

# External Links
$(_doc_external("Sys/PetscDLAddr"))
"""
function PetscDLAddr(petsclib::PetscLibType, func::PetscVoidFn, name::Vector{Cchar}) end

@for_petsc function PetscDLAddr(petsclib::$UnionPetscLib, func::PetscVoidFn, name::Vector{Cchar} )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscDLAddr, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscVoidFn}, Ptr{Ptr{Cchar}}),
               func, name_,
              )


	return nothing
end 

"""
	PetscDemangleSymbol(petsclib::PetscLibType,mangledName::Vector{Cchar}, name::Cchar) 

# External Links
$(_doc_external("Sys/PetscDemangleSymbol"))
"""
function PetscDemangleSymbol(petsclib::PetscLibType, mangledName::Vector{Cchar}, name::Cchar) end

@for_petsc function PetscDemangleSymbol(petsclib::$UnionPetscLib, mangledName::Vector{Cchar}, name::Cchar )

    @chk ccall(
               (:PetscDemangleSymbol, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar),
               mangledName, name,
              )


	return nothing
end 

"""
	PetscErrorPrintfInitialize(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscErrorPrintfInitialize"))
"""
function PetscErrorPrintfInitialize(petsclib::PetscLibType) end

@for_petsc function PetscErrorPrintfInitialize(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscErrorPrintfInitialize, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscFPTrapPush(petsclib::PetscLibType,trap::PetscFPTrap) 
push a floating point trapping mode, restored using `PetscFPTrapPop()`

Not Collective

Input Parameter:
- `trap` - `PETSC_FP_TRAP_ON` or `PETSC_FP_TRAP_OFF` or any of the values passable to `PetscSetFPTrap()`

Level: advanced

-seealso: `PetscFPTrapPop()`, `PetscSetFPTrap()`, `PetscDetermineInitialFPTrap()`

# External Links
$(_doc_external("Sys/PetscFPTrapPush"))
"""
function PetscFPTrapPush(petsclib::PetscLibType, trap::PetscFPTrap) end

@for_petsc function PetscFPTrapPush(petsclib::$UnionPetscLib, trap::PetscFPTrap )

    @chk ccall(
               (:PetscFPTrapPush, $petsc_library),
               PetscErrorCode,
               (PetscFPTrap,),
               trap,
              )


	return nothing
end 

"""
	PetscFPTrapPop(petsclib::PetscLibType) 
push a floating point trapping mode, to be restored using `PetscFPTrapPop()`

Not Collective

Level: advanced

-seealso: `PetscFPTrapPush()`, `PetscSetFPTrap()`, `PetscDetermineInitialFPTrap()`

# External Links
$(_doc_external("Sys/PetscFPTrapPop"))
"""
function PetscFPTrapPop(petsclib::PetscLibType) end

@for_petsc function PetscFPTrapPop(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFPTrapPop, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscSetFPTrap(petsclib::PetscLibType,flag::PetscFPTrap) 
Enables traps/exceptions on common floating point errors. This option may not work on certain systems or only a
subset of exceptions may be trapable.

Not Collective

Input Parameter:
- `flag`  - values are
-seealso: `PetscFPTrapPush()`, `PetscFPTrapPop()`, `PetscDetermineInitialFPTrap()`

# External Links
$(_doc_external("Sys/PetscSetFPTrap"))
"""
function PetscSetFPTrap(petsclib::PetscLibType, flag::PetscFPTrap) end

@for_petsc function PetscSetFPTrap(petsclib::$UnionPetscLib, flag::PetscFPTrap )

    @chk ccall(
               (:PetscSetFPTrap, $petsc_library),
               PetscErrorCode,
               (PetscFPTrap,),
               flag,
              )


	return nothing
end 

"""
	PetscDetermineInitialFPTrap(petsclib::PetscLibType) 
Attempts to determine the floating point trapping that exists when `PetscInitialize()` is called

Not Collective

-seealso: `PetscFPTrapPush()`, `PetscFPTrapPop()`, `PetscDetermineInitialFPTrap()`

# External Links
$(_doc_external("Sys/PetscDetermineInitialFPTrap"))
"""
function PetscDetermineInitialFPTrap(petsclib::PetscLibType) end

@for_petsc function PetscDetermineInitialFPTrap(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDetermineInitialFPTrap, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscCheckPointerSetIntensity(petsclib::PetscLibType,intensity::PetscInt) 
Set the intensity of debug pointer checks

Not Collective

Input Parameter:
- `intensity` - how much to check pointers for validity

Options Database Key:
- `-check_pointer_intensity` - intensity (0, 1, or 2)

Level: advanced

-seealso: `PetscCheckPointer()`, `PetscFunctionBeginHot()`

# External Links
$(_doc_external("Sys/PetscCheckPointerSetIntensity"))
"""
function PetscCheckPointerSetIntensity(petsclib::PetscLibType, intensity::PetscInt) end

@for_petsc function PetscCheckPointerSetIntensity(petsclib::$UnionPetscLib, intensity::$PetscInt )

    @chk ccall(
               (:PetscCheckPointerSetIntensity, $petsc_library),
               PetscErrorCode,
               ($PetscInt,),
               intensity,
              )


	return nothing
end 

"""
	PetscSetDebugTerminal(petsclib::PetscLibType,terminal::Vector{Cchar}) 
Sets the terminal to use for debugging.

Not Collective; No Fortran Support

Input Parameter:
- `terminal` - name of terminal and any flags required to execute a program.
For example "xterm", "urxvt -e", "gnome-terminal -x".
On Apple macOS you can use "Terminal" (note the capital T)

Options Database Key:
- `-debug_terminal terminal` - use this terminal instead of the default

Level: developer

-seealso: `PetscSetDebugger()`, `PetscAttachDebugger()`

# External Links
$(_doc_external("Sys/PetscSetDebugTerminal"))
"""
function PetscSetDebugTerminal(petsclib::PetscLibType, terminal::Vector{Cchar}) end

@for_petsc function PetscSetDebugTerminal(petsclib::$UnionPetscLib, terminal::Vector{Cchar} )

    @chk ccall(
               (:PetscSetDebugTerminal, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               terminal,
              )


	return nothing
end 

"""
	PetscSetDebugger(petsclib::PetscLibType,debugger::Vector{Cchar}, usedebugterminal::PetscBool) 
Sets options associated with the debugger.

Not Collective; No Fortran Support

Input Parameters:
- `debugger`         - name of debugger, which should be in your path,
usually "lldb", "dbx", "gdb", "cuda-gdb", "idb", "xxgdb", "kdgb" or "ddd". Also, HP-UX
supports "xdb", and IBM rs6000 supports "xldb".

- `usedebugterminal` - flag to indicate debugger window, set to either `PETSC_TRUE` (to indicate
debugger should be started in a new terminal window) or `PETSC_FALSE` (to start debugger
in initial window (the option `PETSC_FALSE` makes no sense when using more
than one MPI process.)

Level: developer

-seealso: `PetscAttachDebugger()`, `PetscAttachDebuggerErrorHandler()`, `PetscSetDebugTerminal()`

# External Links
$(_doc_external("Sys/PetscSetDebugger"))
"""
function PetscSetDebugger(petsclib::PetscLibType, debugger::Vector{Cchar}, usedebugterminal::PetscBool) end

@for_petsc function PetscSetDebugger(petsclib::$UnionPetscLib, debugger::Vector{Cchar}, usedebugterminal::PetscBool )

    @chk ccall(
               (:PetscSetDebugger, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscBool),
               debugger, usedebugterminal,
              )


	return nothing
end 

"""
	PetscSetDefaultDebugger(petsclib::PetscLibType) 
Causes PETSc to use its default debugger and output terminal

Not Collective, No Fortran Support

Level: developer

-seealso: `PetscSetDebugger()`, `PetscSetDebuggerFromString()`, `PetscAttachDebugger()`

# External Links
$(_doc_external("Sys/PetscSetDefaultDebugger"))
"""
function PetscSetDefaultDebugger(petsclib::PetscLibType) end

@for_petsc function PetscSetDefaultDebugger(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscSetDefaultDebugger, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscSetDebuggerFromString(petsclib::PetscLibType,string::Vector{Cchar}) 
Set the complete path for the
debugger for PETSc to use.

Not Collective

Input Parameter:
- `string` - the name of the debugger, for example "gdb"

Level: developer

-seealso: `PetscSetDebugger()`, `PetscSetDefaultDebugger()`, `PetscAttachDebugger()`

# External Links
$(_doc_external("Sys/PetscSetDebuggerFromString"))
"""
function PetscSetDebuggerFromString(petsclib::PetscLibType, string::Vector{Cchar}) end

@for_petsc function PetscSetDebuggerFromString(petsclib::$UnionPetscLib, string::Vector{Cchar} )

    @chk ccall(
               (:PetscSetDebuggerFromString, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               string,
              )


	return nothing
end 

"""
	PetscWaitOnError(petsclib::PetscLibType) 
If an error is detected and the process would normally exit the main program with `MPI_Abort()` sleep instead
of exiting.

Not Collective

Level: advanced

-seealso: `PetscSetDebugger()`, `PetscAttachDebugger()`

# External Links
$(_doc_external("Sys/PetscWaitOnError"))
"""
function PetscWaitOnError(petsclib::PetscLibType) end

@for_petsc function PetscWaitOnError(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscWaitOnError, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscAttachDebugger(petsclib::PetscLibType) 
Attaches the debugger to the running process.

Not Collective

Options Database Keys:
- `-start_in_debugger [noxterm,lldb or gdb]` - Set debugger debug_terminal xterm or Terminal (for Apple)
- `-display name`                            - XDisplay to open xterm in
- `-debugger_ranks m,n`                      - Which MPI ranks on which to start the debugger, defaults to all
- `-stop_for_debugger`                       - Print a message on how to attach the process with a debugger and then wait for the user to attach
- `-debugger_pause <secs>`                   - Wait <secs> before attaching the debugger. This is useful for slow connections
that take a long time for the Terminal window or xterm to start up.

Level: advanced

-seealso: `PetscSetDebugger()`, `PetscSetDefaultDebugger()`, `PetscSetDebugTerminal()`, `PetscAttachDebuggerErrorHandler()`, `PetscStopForDebugger()`

# External Links
$(_doc_external("Sys/PetscAttachDebugger"))
"""
function PetscAttachDebugger(petsclib::PetscLibType) end

@for_petsc function PetscAttachDebugger(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscAttachDebugger, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscStopForDebugger(petsclib::PetscLibType) 
Prints a message to the screen indicating how to
attach to the process with the debugger and then waits for the
debugger to attach.

Not Collective, No Fortran Support

Options Database Key:
- `-stop_for_debugger` - will stop for you to attach the debugger when `PetscInitialize()` is called

Level: developer

-seealso: `PetscSetDebugger()`, `PetscAttachDebugger()`

# External Links
$(_doc_external("Sys/PetscStopForDebugger"))
"""
function PetscStopForDebugger(petsclib::PetscLibType) end

@for_petsc function PetscStopForDebugger(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscStopForDebugger, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscPushErrorHandler(petsclib::PetscLibType,handler::external, ctx::Cvoid) 
Sets a routine to be called on detection of errors.

Not Collective, No Fortran Support

Input Parameters:
- `handler` - error handler routine
- `ctx`     - optional handler context that contains information needed by the handler (for
example file pointers for error messages etc.)

Calling sequence of `handler`:
- `comm` - communicator over which error occurred
- `line` - the line number of the error (usually indicated by `__LINE__` in the calling routine)
- `file` - the file in which the error was detected (usually indicated by `__FILE__` in the calling routine)
- `fun`  - the function name of the calling routine
- `n`    - the generic error number (see list defined in include/petscerror.h)
- `p`    - `PETSC_ERROR_INITIAL` if error just detected, otherwise `PETSC_ERROR_REPEAT`
- `mess` - an error text string, usually just printed to the screen
- `ctx`  - the error handler context

Options Database Keys:
- `-on_error_attach_debugger <noxterm,lldb or gdb>` - starts up the debugger if an error occurs
- `-on_error_abort`                                 - aborts the program if an error occurs

Level: intermediate

-seealso: `PetscPopErrorHandler()`, `PetscAttachDebuggerErrorHandler()`, `PetscAbortErrorHandler()`, `PetscTraceBackErrorHandler()`, `PetscPushSignalHandler()`,
`PetscErrorType`, `PETSC_ERROR_INITIAL`, `PETSC_ERROR_REPEAT`, `PetscErrorCode`

# External Links
$(_doc_external("Sys/PetscPushErrorHandler"))
"""
function PetscPushErrorHandler(petsclib::PetscLibType, handler::external, ctx::Cvoid) end

@for_petsc function PetscPushErrorHandler(petsclib::$UnionPetscLib, handler::external, ctx::Cvoid )

    @chk ccall(
               (:PetscPushErrorHandler, $petsc_library),
               PetscErrorCode,
               (external, Ptr{Cvoid}),
               handler, ctx,
              )


	return nothing
end 

"""
	PetscPopErrorHandler(petsclib::PetscLibType) 
Removes the latest error handler that was
pushed with `PetscPushErrorHandler()`.

Not Collective

Level: intermediate

-seealso: `PetscPushErrorHandler()`

# External Links
$(_doc_external("Sys/PetscPopErrorHandler"))
"""
function PetscPopErrorHandler(petsclib::PetscLibType) end

@for_petsc function PetscPopErrorHandler(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPopErrorHandler, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscIntViewNumColumns(petsclib::PetscLibType,N::PetscInt, Ncol::PetscInt, idx::Vector{PetscInt}, viewer::PetscViewer) 
Prints an array of integers; useful for debugging.

Collective

Input Parameters:
- `N`      - number of integers in array
- `Ncol`   - number of integers to print per row
- `idx`    - array of integers
- `viewer` - location to print array, `PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF` or 0

Level: intermediate

-seealso: `PetscViewer`, `PetscIntView()`, `PetscRealView()`

# External Links
$(_doc_external("Sys/PetscIntViewNumColumns"))
"""
function PetscIntViewNumColumns(petsclib::PetscLibType, N::PetscInt, Ncol::PetscInt, idx::Vector{PetscInt}, viewer::PetscViewer) end

@for_petsc function PetscIntViewNumColumns(petsclib::$UnionPetscLib, N::$PetscInt, Ncol::$PetscInt, idx::Vector{$PetscInt}, viewer::PetscViewer )

    @chk ccall(
               (:PetscIntViewNumColumns, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}, PetscViewer),
               N, Ncol, idx, viewer,
              )


	return nothing
end 

"""
	PetscRealViewNumColumns(petsclib::PetscLibType,N::PetscInt, Ncol::PetscInt, idx::Vector{PetscReal}, viewer::PetscViewer) 
Prints an array of doubles; useful for debugging.

Collective

Input Parameters:
- `N`      - number of `PetscReal` in array
- `Ncol`   - number of `PetscReal` to print per row
- `idx`    - array of `PetscReal`
- `viewer` - location to print array, `PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF` or 0

Level: intermediate

-seealso: `PetscViewer`, `PetscRealView()`, `PetscIntView()`

# External Links
$(_doc_external("Sys/PetscRealViewNumColumns"))
"""
function PetscRealViewNumColumns(petsclib::PetscLibType, N::PetscInt, Ncol::PetscInt, idx::Vector{PetscReal}, viewer::PetscViewer) end

@for_petsc function PetscRealViewNumColumns(petsclib::$UnionPetscLib, N::$PetscInt, Ncol::$PetscInt, idx::Vector{$PetscReal}, viewer::PetscViewer )

    @chk ccall(
               (:PetscRealViewNumColumns, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, PetscViewer),
               N, Ncol, idx, viewer,
              )


	return nothing
end 

"""
	PetscScalarViewNumColumns(petsclib::PetscLibType,N::PetscInt, Ncol::PetscInt, idx::Vector{PetscScalar}, viewer::PetscViewer) 
Prints an array of doubles; useful for debugging.

Collective

Input Parameters:
- `N`      - number of `PetscScalar` in array
- `Ncol`   - number of `PetscScalar` to print per row
- `idx`    - array of `PetscScalar`
- `viewer` - location to print array, `PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF` or 0

Level: intermediate

-seealso: `PetscViewer`, `PetscRealView()`, `PetscScalarView()`, `PetscIntView()`

# External Links
$(_doc_external("Sys/PetscScalarViewNumColumns"))
"""
function PetscScalarViewNumColumns(petsclib::PetscLibType, N::PetscInt, Ncol::PetscInt, idx::Vector{PetscScalar}, viewer::PetscViewer) end

@for_petsc function PetscScalarViewNumColumns(petsclib::$UnionPetscLib, N::$PetscInt, Ncol::$PetscInt, idx::Vector{$PetscScalar}, viewer::PetscViewer )

    @chk ccall(
               (:PetscScalarViewNumColumns, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscScalar}, PetscViewer),
               N, Ncol, idx, viewer,
              )


	return nothing
end 

"""
	PetscIntView(petsclib::PetscLibType,N::PetscInt, idx::Vector{PetscInt}, viewer::PetscViewer) 
Prints an array of integers; useful for debugging.

Collective

Input Parameters:
- `N`      - number of integers in array
- `idx`    - array of integers
- `viewer` - location to print array, `PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF` or 0

Level: intermediate

-seealso: `PetscViewer`, `PetscIntViewNumColumns()`, `PetscRealView()`

# External Links
$(_doc_external("Sys/PetscIntView"))
"""
function PetscIntView(petsclib::PetscLibType, N::PetscInt, idx::Vector{PetscInt}, viewer::PetscViewer) end

@for_petsc function PetscIntView(petsclib::$UnionPetscLib, N::$PetscInt, idx::Vector{$PetscInt}, viewer::PetscViewer )

    @chk ccall(
               (:PetscIntView, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, PetscViewer),
               N, idx, viewer,
              )


	return nothing
end 

"""
	PetscRealView(petsclib::PetscLibType,N::PetscInt, idx::Vector{PetscReal}, viewer::PetscViewer) 
Prints an array of doubles; useful for debugging.

Collective

Input Parameters:
- `N`      - number of `PetscReal` in array
- `idx`    - array of `PetscReal`
- `viewer` - location to print array, `PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF` or 0

Level: intermediate

-seealso: `PetscViewer`, `PetscIntView()`

# External Links
$(_doc_external("Sys/PetscRealView"))
"""
function PetscRealView(petsclib::PetscLibType, N::PetscInt, idx::Vector{PetscReal}, viewer::PetscViewer) end

@for_petsc function PetscRealView(petsclib::$UnionPetscLib, N::$PetscInt, idx::Vector{$PetscReal}, viewer::PetscViewer )

    @chk ccall(
               (:PetscRealView, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, PetscViewer),
               N, idx, viewer,
              )


	return nothing
end 

"""
	PetscScalarView(petsclib::PetscLibType,N::PetscInt, idx::Vector{PetscScalar}, viewer::PetscViewer) 
Prints an array of `PetscScalar`; useful for debugging.

Collective

Input Parameters:
- `N`      - number of scalars in array
- `idx`    - array of scalars
- `viewer` - location to print array, `PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF` or 0

Level: intermediate

-seealso: `PetscViewer`, `PetscIntView()`, `PetscRealView()`

# External Links
$(_doc_external("Sys/PetscScalarView"))
"""
function PetscScalarView(petsclib::PetscLibType, N::PetscInt, idx::Vector{PetscScalar}, viewer::PetscViewer) end

@for_petsc function PetscScalarView(petsclib::$UnionPetscLib, N::$PetscInt, idx::Vector{$PetscScalar}, viewer::PetscViewer )

    @chk ccall(
               (:PetscScalarView, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscScalar}, PetscViewer),
               N, idx, viewer,
              )


	return nothing
end 

"""
	PetscSignalHandlerDefault(petsclib::PetscLibType,sig::Cint, ptr::Cvoid) 
Default signal handler.

Not Collective

Input Parameters:
- `sig` - signal value
- `ptr` - unused pointer

Level: advanced

-seealso: [](sec_errors), `PetscPushSignalHandler()`

# External Links
$(_doc_external("Sys/PetscSignalHandlerDefault"))
"""
function PetscSignalHandlerDefault(petsclib::PetscLibType, sig::Cint, ptr::Cvoid) end

@for_petsc function PetscSignalHandlerDefault(petsclib::$UnionPetscLib, sig::Cint, ptr::Cvoid )

    @chk ccall(
               (:PetscSignalHandlerDefault, $petsc_library),
               PetscErrorCode,
               (Cint, Ptr{Cvoid}),
               sig, ptr,
              )


	return nothing
end 

"""
	PetscPushSignalHandler(petsclib::PetscLibType,routine::external, ctx::Cvoid) 
Catches the usual fatal errors and
calls a user-provided routine.

Not Collective, No Fortran Support

Input Parameters:
- `routine` - routine to call when a signal is received
- `ctx`     - optional context needed by the routine

Level: developer

-seealso: [](sec_errors), `PetscPopSignalHandler()`, `PetscSignalHandlerDefault()`, `PetscPushErrorHandler()`

# External Links
$(_doc_external("Sys/PetscPushSignalHandler"))
"""
function PetscPushSignalHandler(petsclib::PetscLibType, routine::external, ctx::Cvoid) end

@for_petsc function PetscPushSignalHandler(petsclib::$UnionPetscLib, routine::external, ctx::Cvoid )

    @chk ccall(
               (:PetscPushSignalHandler, $petsc_library),
               PetscErrorCode,
               (external, Ptr{Cvoid}),
               routine, ctx,
              )


	return nothing
end 

"""
	PetscPopSignalHandler(petsclib::PetscLibType) 
Removes the last signal handler that was pushed.
If no signal handlers are left on the stack it will remove the PETSc signal handler.
(That is PETSc will no longer catch signals).

Not Collective

Level: developer

-seealso: [](sec_errors), `PetscPushSignalHandler()`

# External Links
$(_doc_external("Sys/PetscPopSignalHandler"))
"""
function PetscPopSignalHandler(petsclib::PetscLibType) end

@for_petsc function PetscPopSignalHandler(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscPopSignalHandler, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscAddLogDouble(petsclib::PetscLibType,tot::PetscLogDouble, tot_th::PetscLogDouble, tmp::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscAddLogDouble"))
"""
function PetscAddLogDouble(petsclib::PetscLibType, tot::PetscLogDouble, tot_th::PetscLogDouble, tmp::PetscLogDouble) end

@for_petsc function PetscAddLogDouble(petsclib::$UnionPetscLib, tot::PetscLogDouble, tot_th::PetscLogDouble, tmp::PetscLogDouble )

    @chk ccall(
               (:PetscAddLogDouble, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble}, Ptr{PetscLogDouble}, PetscLogDouble),
               tot, tot_th, tmp,
              )


	return nothing
end 

"""
	PetscAddLogDoubleCnt(petsclib::PetscLibType,cnt::PetscLogDouble, tot::PetscLogDouble, cnt_th::PetscLogDouble, tot_th::PetscLogDouble, tmp::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscAddLogDoubleCnt"))
"""
function PetscAddLogDoubleCnt(petsclib::PetscLibType, cnt::PetscLogDouble, tot::PetscLogDouble, cnt_th::PetscLogDouble, tot_th::PetscLogDouble, tmp::PetscLogDouble) end

@for_petsc function PetscAddLogDoubleCnt(petsclib::$UnionPetscLib, cnt::PetscLogDouble, tot::PetscLogDouble, cnt_th::PetscLogDouble, tot_th::PetscLogDouble, tmp::PetscLogDouble )

    @chk ccall(
               (:PetscAddLogDoubleCnt, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble}, Ptr{PetscLogDouble}, Ptr{PetscLogDouble}, Ptr{PetscLogDouble}, PetscLogDouble),
               cnt, tot, cnt_th, tot_th, tmp,
              )


	return nothing
end 

"""
	PetscLogGetDefaultHandler(petsclib::PetscLibType,handler::PetscLogHandler) 
Get the default log handler if it is running.

Not collective

Output Parameter:
- `handler` - the default `PetscLogHandler`, or `NULL` if it is not running.

Level: developer

-seealso: [](ch_profiling)

# External Links
$(_doc_external("Sys/PetscLogGetDefaultHandler"))
"""
function PetscLogGetDefaultHandler(petsclib::PetscLibType, handler::PetscLogHandler) end

@for_petsc function PetscLogGetDefaultHandler(petsclib::$UnionPetscLib, handler::PetscLogHandler )

    @chk ccall(
               (:PetscLogGetDefaultHandler, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogHandler},),
               handler,
              )


	return nothing
end 

"""
	PetscLogGetState(petsclib::PetscLibType,state::PetscLogState) 
Get the `PetscLogState` for PETSc's global logging, used
by all default log handlers (`PetscLogDefaultBegin()`,
`PetscLogNestedBegin()`, `PetscLogTraceBegin()`, `PetscLogMPEBegin()`,
`PetscLogPerfstubsBegin()`).

Collective on `PETSC_COMM_WORLD`

Output Parameter:
- `state` - The `PetscLogState` changed by registrations (such as
`PetscLogEventRegister()`) and actions (such as `PetscLogEventBegin()` or
`PetscLogStagePush()`), or `NULL` if logging is not active

Level: developer

-seealso: [](ch_profiling), `PetscLogState`

# External Links
$(_doc_external("Sys/PetscLogGetState"))
"""
function PetscLogGetState(petsclib::PetscLibType, state::PetscLogState) end

@for_petsc function PetscLogGetState(petsclib::$UnionPetscLib, state::PetscLogState )

    @chk ccall(
               (:PetscLogGetState, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogState},),
               state,
              )


	return nothing
end 

"""
	isActive::PetscBool = PetscLogIsActive(petsclib::PetscLibType) 
Check if logging (profiling) is currently in progress.

Not Collective

Output Parameter:
- `isActive` - `PETSC_TRUE` if logging is in progress, `PETSC_FALSE` otherwise

Level: beginner

-seealso: [](ch_profiling), `PetscLogDefaultBegin()`

# External Links
$(_doc_external("Sys/PetscLogIsActive"))
"""
function PetscLogIsActive(petsclib::PetscLibType) end

@for_petsc function PetscLogIsActive(petsclib::$UnionPetscLib)
	isActive_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLogIsActive, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               isActive_,
              )

	isActive = isActive_[]

	return isActive
end 

"""
	PetscLogDefaultBegin(petsclib::PetscLibType) 
Turns on logging (profiling) of PETSc code using the default log handler (profiler). This logs time, flop
rates, and object creation and should not slow programs down too much.

Logically Collective on `PETSC_COMM_WORLD`

Options Database Key:
- `-log_view [viewertype:filename:viewerformat]` - Prints summary of flop and timing (profiling) information to the
screen (for PETSc configured with `--with-log=1` (which is the default)).
This option must be provided before `PetscInitialize()`.

-seealso: [](ch_profiling), `PetscLogDump()`, `PetscLogView()`, `PetscLogTraceBegin()`

# External Links
$(_doc_external("Sys/PetscLogDefaultBegin"))
"""
function PetscLogDefaultBegin(petsclib::PetscLibType) end

@for_petsc function PetscLogDefaultBegin(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogDefaultBegin, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogTraceBegin(petsclib::PetscLibType,file::Libc.FILE) 
Begins trace logging.  Every time a PETSc event
begins or ends, the event name is printed.

Logically Collective on `PETSC_COMM_WORLD`, No Fortran Support

Input Parameter:
- `file` - The file to print trace in (e.g. stdout)

Options Database Key:
- `-log_trace [filename]` - Begins `PetscLogTraceBegin()`

Level: intermediate

-seealso: [](ch_profiling), `PetscLogDump()`, `PetscLogView()`, `PetscLogDefaultBegin()`

# External Links
$(_doc_external("Sys/PetscLogTraceBegin"))
"""
function PetscLogTraceBegin(petsclib::PetscLibType, file::Libc.FILE) end

@for_petsc function PetscLogTraceBegin(petsclib::$UnionPetscLib, file::Libc.FILE )

    @chk ccall(
               (:PetscLogTraceBegin, $petsc_library),
               PetscErrorCode,
               (Ptr{Libc.FILE},),
               file,
              )


	return nothing
end 

"""
	PetscLogNestedBegin(petsclib::PetscLibType) 
Turns on nested logging of objects and events. This logs flop
rates and object creation and should not slow programs down too much.

Logically Collective on `PETSC_COMM_WORLD`, No Fortran Support

Options Database Keys:
- `-log_view :filename.xml:ascii_xml` - Prints an XML summary of flop and timing information to the file

-seealso: `PetscLogDump()`, `PetscLogView()`, `PetscLogTraceBegin()`, `PetscLogDefaultBegin()`

# External Links
$(_doc_external("Sys/PetscLogNestedBegin"))
"""
function PetscLogNestedBegin(petsclib::PetscLibType) end

@for_petsc function PetscLogNestedBegin(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogNestedBegin, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogLegacyCallbacksBegin(petsclib::PetscLibType,PetscLogPLB::external, PetscLogPLE::external, PetscLogPHC::external, PetscLogPHD::external) 
Create and start a log handler from callbacks
matching the now deprecated function pointers `PetscLogPLB`, `PetscLogPLE`,
`PetscLogPHC`, `PetscLogPHD`.

Logically Collective on `PETSC_COMM_WORLD`

Input Parameters:
- `PetscLogPLB` - A callback that will be executed by `PetscLogEventBegin()` (or `NULL`)
- `PetscLogPLE` - A callback that will be executed by `PetscLogEventEnd()` (or `NULL`)
- `PetscLogPHC` - A callback that will be executed by `PetscLogObjectCreate()` (or `NULL`)
- `PetscLogPHD` - A callback that will be executed by `PetscLogObjectCreate()` (or `NULL`)

Calling sequence of `PetscLogPLB`:
- `e`  - a `PetscLogEvent` that is beginning
- `_i` - deprecated, unused
- `o1` - a `PetscObject` associated with `e` (or `NULL`)
- `o2` - a `PetscObject` associated with `e` (or `NULL`)
- `o3` - a `PetscObject` associated with `e` (or `NULL`)
- `o4` - a `PetscObject` associated with `e` (or `NULL`)

Calling sequence of `PetscLogPLE`:
- `e`  - a `PetscLogEvent` that is beginning
- `_i` - deprecated, unused
- `o1` - a `PetscObject` associated with `e` (or `NULL`)
- `o2` - a `PetscObject` associated with `e` (or `NULL`)
- `o3` - a `PetscObject` associated with `e` (or `NULL`)
- `o4` - a `PetscObject` associated with `e` (or `NULL`)

Calling sequence of `PetscLogPHC`:
- `o` - a `PetscObject` that has just been created

Calling sequence of `PetscLogPHD`:
- `o` - a `PetscObject` that is about to be destroyed

Level: advanced

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerStart()`, `PetscLogState`

# External Links
$(_doc_external("Sys/PetscLogLegacyCallbacksBegin"))
"""
function PetscLogLegacyCallbacksBegin(petsclib::PetscLibType, PetscLogPLB::external, PetscLogPLE::external, PetscLogPHC::external, PetscLogPHD::external) end

@for_petsc function PetscLogLegacyCallbacksBegin(petsclib::$UnionPetscLib, PetscLogPLB::external, PetscLogPLE::external, PetscLogPHC::external, PetscLogPHD::external )

    @chk ccall(
               (:PetscLogLegacyCallbacksBegin, $petsc_library),
               PetscErrorCode,
               (external, external, external, external),
               PetscLogPLB, PetscLogPLE, PetscLogPHC, PetscLogPHD,
              )


	return nothing
end 

"""
	PetscLogMPEBegin(petsclib::PetscLibType) 
Turns on MPE logging of events. This creates large log files and slows the
program down.

Collective on `PETSC_COMM_WORLD`, No Fortran Support

Options Database Key:
- `-log_mpe` - Prints extensive log information

Level: advanced

-seealso: [](ch_profiling), `PetscLogDump()`, `PetscLogDefaultBegin()`, `PetscLogEventActivate()`,
`PetscLogEventDeactivate()`

# External Links
$(_doc_external("Sys/PetscLogMPEBegin"))
"""
function PetscLogMPEBegin(petsclib::PetscLibType) end

@for_petsc function PetscLogMPEBegin(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogMPEBegin, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogPerfstubsBegin(petsclib::PetscLibType) 
Turns on logging of events using the perfstubs interface.

Collective on `PETSC_COMM_WORLD`, No Fortran Support

Options Database Key:
- `-log_perfstubs` - use an external log handler through the perfstubs interface

Level: advanced

-seealso: [](ch_profiling), `PetscLogDefaultBegin()`, `PetscLogEventActivate()`

# External Links
$(_doc_external("Sys/PetscLogPerfstubsBegin"))
"""
function PetscLogPerfstubsBegin(petsclib::PetscLibType) end

@for_petsc function PetscLogPerfstubsBegin(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogPerfstubsBegin, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogActions(petsclib::PetscLibType,flag::PetscBool) 
Determines whether actions are logged for the default log handler.

Not Collective

Input Parameter:
- `flag` - `PETSC_TRUE` if actions are to be logged

Options Database Key:
- `-log_exclude_actions` - (deprecated) Does nothing
- `-log_include_actions` - Turn on action logging

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogActions"))
"""
function PetscLogActions(petsclib::PetscLibType, flag::PetscBool) end

@for_petsc function PetscLogActions(petsclib::$UnionPetscLib, flag::PetscBool )

    @chk ccall(
               (:PetscLogActions, $petsc_library),
               PetscErrorCode,
               (PetscBool,),
               flag,
              )


	return nothing
end 

"""
	PetscLogObjects(petsclib::PetscLibType,flag::PetscBool) 
Determines whether objects are logged for the graphical viewer.

Not Collective

Input Parameter:
- `flag` - `PETSC_TRUE` if objects are to be logged

Options Database Key:
- `-log_exclude_objects` - (deprecated) Does nothing
- `-log_include_objects` - Turns on object logging

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogObjects"))
"""
function PetscLogObjects(petsclib::PetscLibType, flag::PetscBool) end

@for_petsc function PetscLogObjects(petsclib::$UnionPetscLib, flag::PetscBool )

    @chk ccall(
               (:PetscLogObjects, $petsc_library),
               PetscErrorCode,
               (PetscBool,),
               flag,
              )


	return nothing
end 

"""
	PetscLogStageRegister(petsclib::PetscLibType,sname::Vector{Cchar}, stage::PetscLogStage) 
Attaches a character string name to a logging stage.

Not Collective

Input Parameter:
- `sname` - The name to associate with that stage

Output Parameter:
- `stage` - The stage number or -1 if logging is not active (`PetscLogIsActive()`).

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStagePush()`, `PetscLogStagePop()`

# External Links
$(_doc_external("Sys/PetscLogStageRegister"))
"""
function PetscLogStageRegister(petsclib::PetscLibType, sname::Vector{Cchar}, stage::PetscLogStage) end

@for_petsc function PetscLogStageRegister(petsclib::$UnionPetscLib, sname::Vector{Cchar}, stage::PetscLogStage )

    @chk ccall(
               (:PetscLogStageRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscLogStage}),
               sname, stage,
              )


	return nothing
end 

"""
	PetscLogStagePush(petsclib::PetscLibType,stage::PetscLogStage) 
This function pushes a stage on the logging stack. Events started and stopped until `PetscLogStagePop()` will be associated with the stage

Not Collective

Input Parameter:
- `stage` - The stage on which to log

-seealso: [](ch_profiling), `PetscLogStagePop()`, `PetscLogStageRegister()`, `PetscBarrier()`

# External Links
$(_doc_external("Sys/PetscLogStagePush"))
"""
function PetscLogStagePush(petsclib::PetscLibType, stage::PetscLogStage) end

@for_petsc function PetscLogStagePush(petsclib::$UnionPetscLib, stage::PetscLogStage )

    @chk ccall(
               (:PetscLogStagePush, $petsc_library),
               PetscErrorCode,
               (PetscLogStage,),
               stage,
              )


	return nothing
end 

"""
	PetscLogStagePop(petsclib::PetscLibType) 
This function pops a stage from the logging stack that was pushed with `PetscLogStagePush()`

Not Collective

-seealso: [](ch_profiling), `PetscLogStagePush()`, `PetscLogStageRegister()`, `PetscBarrier()`

# External Links
$(_doc_external("Sys/PetscLogStagePop"))
"""
function PetscLogStagePop(petsclib::PetscLibType) end

@for_petsc function PetscLogStagePop(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogStagePop, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogStageSetActive(petsclib::PetscLibType,stage::PetscLogStage, isActive::PetscBool) 
Sets if a stage is used for `PetscLogEventBegin()` and `PetscLogEventEnd()`.

Not Collective

Input Parameters:
- `stage`    - The stage
- `isActive` - The activity flag, `PETSC_TRUE` for logging, else `PETSC_FALSE` (defaults to `PETSC_TRUE`)

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStageRegister()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscPreLoadBegin()`, `PetscPreLoadEnd()`, `PetscPreLoadStage()`

# External Links
$(_doc_external("Sys/PetscLogStageSetActive"))
"""
function PetscLogStageSetActive(petsclib::PetscLibType, stage::PetscLogStage, isActive::PetscBool) end

@for_petsc function PetscLogStageSetActive(petsclib::$UnionPetscLib, stage::PetscLogStage, isActive::PetscBool )

    @chk ccall(
               (:PetscLogStageSetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, PetscBool),
               stage, isActive,
              )


	return nothing
end 

"""
	isActive::PetscBool = PetscLogStageGetActive(petsclib::PetscLibType,stage::PetscLogStage) 
Checks if a stage is used for `PetscLogEventBegin()` and `PetscLogEventEnd()`.

Not Collective

Input Parameter:
- `stage` - The stage

Output Parameter:
- `isActive` - The activity flag, `PETSC_TRUE` for logging, else `PETSC_FALSE` (defaults to `PETSC_TRUE`)

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStageRegister()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscPreLoadBegin()`, `PetscPreLoadEnd()`, `PetscPreLoadStage()`

# External Links
$(_doc_external("Sys/PetscLogStageGetActive"))
"""
function PetscLogStageGetActive(petsclib::PetscLibType, stage::PetscLogStage) end

@for_petsc function PetscLogStageGetActive(petsclib::$UnionPetscLib, stage::PetscLogStage )
	isActive_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLogStageGetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, Ptr{PetscBool}),
               stage, isActive_,
              )

	isActive = isActive_[]

	return isActive
end 

"""
	PetscLogStageSetVisible(petsclib::PetscLibType,stage::PetscLogStage, isVisible::PetscBool) 
Determines stage visibility in `PetscLogView()`

Not Collective

Input Parameters:
- `stage`     - The stage
- `isVisible` - The visibility flag, `PETSC_TRUE` to print, else `PETSC_FALSE` (defaults to `PETSC_TRUE`)

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStageGetVisible()`, `PetscLogStageRegister()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogView()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogStageSetVisible"))
"""
function PetscLogStageSetVisible(petsclib::PetscLibType, stage::PetscLogStage, isVisible::PetscBool) end

@for_petsc function PetscLogStageSetVisible(petsclib::$UnionPetscLib, stage::PetscLogStage, isVisible::PetscBool )

    @chk ccall(
               (:PetscLogStageSetVisible, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, PetscBool),
               stage, isVisible,
              )


	return nothing
end 

"""
	isVisible::PetscBool = PetscLogStageGetVisible(petsclib::PetscLibType,stage::PetscLogStage) 
Returns stage visibility in `PetscLogView()`

Not Collective

Input Parameter:
- `stage` - The stage

Output Parameter:
- `isVisible` - The visibility flag, `PETSC_TRUE` to print, else `PETSC_FALSE` (defaults to `PETSC_TRUE`)

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStageSetVisible()`, `PetscLogStageRegister()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogView()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogStageGetVisible"))
"""
function PetscLogStageGetVisible(petsclib::PetscLibType, stage::PetscLogStage) end

@for_petsc function PetscLogStageGetVisible(petsclib::$UnionPetscLib, stage::PetscLogStage )
	isVisible_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLogStageGetVisible, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, Ptr{PetscBool}),
               stage, isVisible_,
              )

	isVisible = isVisible_[]

	return isVisible
end 

"""
	PetscLogStageGetId(petsclib::PetscLibType,name::Vector{Cchar}, stage::PetscLogStage) 
Returns the stage id when given the stage name.

Not Collective

Input Parameter:
- `name` - The stage name

Output Parameter:
- `stage` - The stage, , or -1 if no stage with that name exists

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStageRegister()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscPreLoadBegin()`, `PetscPreLoadEnd()`, `PetscPreLoadStage()`

# External Links
$(_doc_external("Sys/PetscLogStageGetId"))
"""
function PetscLogStageGetId(petsclib::PetscLibType, name::Vector{Cchar}, stage::PetscLogStage) end

@for_petsc function PetscLogStageGetId(petsclib::$UnionPetscLib, name::Vector{Cchar}, stage::PetscLogStage )

    @chk ccall(
               (:PetscLogStageGetId, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscLogStage}),
               name, stage,
              )


	return nothing
end 

"""
	PetscLogStageGetName(petsclib::PetscLibType,stage::PetscLogStage, name::Vector{Cchar}) 
Returns the stage name when given the stage id.

Not Collective

Input Parameter:
- `stage` - The stage

Output Parameter:
- `name` - The stage name

Level: intermediate

-seealso: [](ch_profiling), `PetscLogStageRegister()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscPreLoadBegin()`, `PetscPreLoadEnd()`, `PetscPreLoadStage()`

# External Links
$(_doc_external("Sys/PetscLogStageGetName"))
"""
function PetscLogStageGetName(petsclib::PetscLibType, stage::PetscLogStage, name::Vector{Cchar}) end

@for_petsc function PetscLogStageGetName(petsclib::$UnionPetscLib, stage::PetscLogStage, name::Vector{Cchar} )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscLogStageGetName, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, Ptr{Ptr{Cchar}}),
               stage, name_,
              )


	return nothing
end 

"""
	PetscLogEventRegister(petsclib::PetscLibType,name::Vector{Cchar}, classid::PetscClassId, event::PetscLogEvent) 
Registers an event name for logging operations

Not Collective

Input Parameters:
- `name`    - The name associated with the event
- `classid` - The classid associated to the class for this event, obtain either with
`PetscClassIdRegister()` or use a predefined one such as `KSP_CLASSID`, `SNES_CLASSID`, the predefined ones
are only available in C code

Output Parameter:
- `event` - The event id for use with `PetscLogEventBegin()` and `PetscLogEventEnd()`.

-seealso: [](ch_profiling), `PetscLogStageRegister()`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogFlops()`,
`PetscLogEventActivate()`, `PetscLogEventDeactivate()`, `PetscClassIdRegister()`

# External Links
$(_doc_external("Sys/PetscLogEventRegister"))
"""
function PetscLogEventRegister(petsclib::PetscLibType, name::Vector{Cchar}, classid::PetscClassId, event::PetscLogEvent) end

@for_petsc function PetscLogEventRegister(petsclib::$UnionPetscLib, name::Vector{Cchar}, classid::PetscClassId, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogEventRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscClassId, Ptr{PetscLogEvent}),
               name, classid, event,
              )


	return nothing
end 

"""
	PetscLogEventSetCollective(petsclib::PetscLibType,event::PetscLogEvent, collective::PetscBool) 
Indicates that a particular event is collective.

Logically Collective

Input Parameters:
- `event`      - The event id
- `collective` - `PetscBool` indicating whether a particular event is collective

Level: developer

-seealso: [](ch_profiling), `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogEventRegister()`

# External Links
$(_doc_external("Sys/PetscLogEventSetCollective"))
"""
function PetscLogEventSetCollective(petsclib::PetscLibType, event::PetscLogEvent, collective::PetscBool) end

@for_petsc function PetscLogEventSetCollective(petsclib::$UnionPetscLib, event::PetscLogEvent, collective::PetscBool )

    @chk ccall(
               (:PetscLogEventSetCollective, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent, PetscBool),
               event, collective,
              )


	return nothing
end 

"""
	PetscLogEventIncludeClass(petsclib::PetscLibType,classid::PetscClassId) 
Activates event logging for a PETSc object class in every stage.

Not Collective

Input Parameter:
- `classid` - The object class, for example `MAT_CLASSID`, `SNES_CLASSID`, etc.

Level: developer

-seealso: [](ch_profiling), `PetscLogEventActivateClass()`, `PetscLogEventDeactivateClass()`, `PetscLogEventActivate()`, `PetscLogEventDeactivate()`

# External Links
$(_doc_external("Sys/PetscLogEventIncludeClass"))
"""
function PetscLogEventIncludeClass(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscLogEventIncludeClass(petsclib::$UnionPetscLib, classid::PetscClassId )

    @chk ccall(
               (:PetscLogEventIncludeClass, $petsc_library),
               PetscErrorCode,
               (PetscClassId,),
               classid,
              )


	return nothing
end 

"""
	PetscLogEventExcludeClass(petsclib::PetscLibType,classid::PetscClassId) 
Deactivates event logging for a PETSc object class in every stage.

Not Collective

Input Parameter:
- `classid` - The object class, for example `MAT_CLASSID`, `SNES_CLASSID`, etc.

Level: developer

-seealso: [](ch_profiling), `PetscLogEventDeactivateClass()`, `PetscLogEventActivateClass()`, `PetscLogEventDeactivate()`, `PetscLogEventActivate()`

# External Links
$(_doc_external("Sys/PetscLogEventExcludeClass"))
"""
function PetscLogEventExcludeClass(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscLogEventExcludeClass(petsclib::$UnionPetscLib, classid::PetscClassId )

    @chk ccall(
               (:PetscLogEventExcludeClass, $petsc_library),
               PetscErrorCode,
               (PetscClassId,),
               classid,
              )


	return nothing
end 

"""
	PetscLogEventActivate(petsclib::PetscLibType,event::PetscLogEvent) 
Indicates that a particular event should be logged.

Not Collective

Input Parameter:
- `event` - The event id

-seealso: [](ch_profiling), `PetscLogEventDeactivate()`, `PetscLogEventDeactivatePush()`, `PetscLogEventDeactivatePop()`

# External Links
$(_doc_external("Sys/PetscLogEventActivate"))
"""
function PetscLogEventActivate(petsclib::PetscLibType, event::PetscLogEvent) end

@for_petsc function PetscLogEventActivate(petsclib::$UnionPetscLib, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogEventActivate, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent,),
               event,
              )


	return nothing
end 

"""
	PetscLogEventDeactivate(petsclib::PetscLibType,event::PetscLogEvent) 
Indicates that a particular event should not be logged.

Not Collective

Input Parameter:
- `event` - The event id

-seealso: [](ch_profiling), `PetscLogEventActivate()`, `PetscLogEventDeactivatePush()`, `PetscLogEventDeactivatePop()`

# External Links
$(_doc_external("Sys/PetscLogEventDeactivate"))
"""
function PetscLogEventDeactivate(petsclib::PetscLibType, event::PetscLogEvent) end

@for_petsc function PetscLogEventDeactivate(petsclib::$UnionPetscLib, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogEventDeactivate, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent,),
               event,
              )


	return nothing
end 

"""
	PetscLogEventDeactivatePush(petsclib::PetscLibType,event::PetscLogEvent) 
Indicates that a particular event should not be logged until `PetscLogEventDeactivatePop()` is called

Not Collective

Input Parameter:
- `event` - The event id

-seealso: [](ch_profiling), `PetscLogEventActivate()`, `PetscLogEventDeactivate()`, `PetscLogEventDeactivatePop()`

# External Links
$(_doc_external("Sys/PetscLogEventDeactivatePush"))
"""
function PetscLogEventDeactivatePush(petsclib::PetscLibType, event::PetscLogEvent) end

@for_petsc function PetscLogEventDeactivatePush(petsclib::$UnionPetscLib, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogEventDeactivatePush, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent,),
               event,
              )


	return nothing
end 

"""
	PetscLogEventDeactivatePop(petsclib::PetscLibType,event::PetscLogEvent) 
Indicates that a particular event should again be logged after the logging was turned off with `PetscLogEventDeactivatePush()`

Not Collective

Input Parameter:
- `event` - The event id

-seealso: [](ch_profiling), `PetscLogEventActivate()`, `PetscLogEventDeactivatePush()`

# External Links
$(_doc_external("Sys/PetscLogEventDeactivatePop"))
"""
function PetscLogEventDeactivatePop(petsclib::PetscLibType, event::PetscLogEvent) end

@for_petsc function PetscLogEventDeactivatePop(petsclib::$UnionPetscLib, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogEventDeactivatePop, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent,),
               event,
              )


	return nothing
end 

"""
	PetscLogEventSetActiveAll(petsclib::PetscLibType,event::PetscLogEvent, isActive::PetscBool) 
Turns on logging of all events

Not Collective

Input Parameters:
- `event`    - The event id
- `isActive` - The activity flag determining whether the event is logged

Level: advanced

-seealso: [](ch_profiling), `PetscLogEventActivate()`, `PetscLogEventDeactivate()`

# External Links
$(_doc_external("Sys/PetscLogEventSetActiveAll"))
"""
function PetscLogEventSetActiveAll(petsclib::PetscLibType, event::PetscLogEvent, isActive::PetscBool) end

@for_petsc function PetscLogEventSetActiveAll(petsclib::$UnionPetscLib, event::PetscLogEvent, isActive::PetscBool )

    @chk ccall(
               (:PetscLogEventSetActiveAll, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent, PetscBool),
               event, isActive,
              )


	return nothing
end 

"""
	PetscLogEventActivateClass(petsclib::PetscLibType,classid::PetscClassId) 
Activates event logging for a PETSc object class for the current stage

Not Collective

Input Parameter:
- `classid` - The event class, for example `MAT_CLASSID`, `SNES_CLASSID`, etc.

Level: developer

-seealso: [](ch_profiling), `PetscLogEventIncludeClass()`, `PetscLogEventExcludeClass()`, `PetscLogEventDeactivateClass()`, `PetscLogEventActivate()`, `PetscLogEventDeactivate()`

# External Links
$(_doc_external("Sys/PetscLogEventActivateClass"))
"""
function PetscLogEventActivateClass(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscLogEventActivateClass(petsclib::$UnionPetscLib, classid::PetscClassId )

    @chk ccall(
               (:PetscLogEventActivateClass, $petsc_library),
               PetscErrorCode,
               (PetscClassId,),
               classid,
              )


	return nothing
end 

"""
	PetscLogEventDeactivateClass(petsclib::PetscLibType,classid::PetscClassId) 
Deactivates event logging for a PETSc object class for the current stage

Not Collective

Input Parameter:
- `classid` - The event class, for example `MAT_CLASSID`, `SNES_CLASSID`, etc.

Level: developer

-seealso: [](ch_profiling), `PetscLogEventIncludeClass()`, `PetscLogEventExcludeClass()`, `PetscLogEventActivateClass()`, `PetscLogEventActivate()`, `PetscLogEventDeactivate()`

# External Links
$(_doc_external("Sys/PetscLogEventDeactivateClass"))
"""
function PetscLogEventDeactivateClass(petsclib::PetscLibType, classid::PetscClassId) end

@for_petsc function PetscLogEventDeactivateClass(petsclib::$UnionPetscLib, classid::PetscClassId )

    @chk ccall(
               (:PetscLogEventDeactivateClass, $petsc_library),
               PetscErrorCode,
               (PetscClassId,),
               classid,
              )


	return nothing
end 

"""
	PetscLogEventSync(petsclib::PetscLibType,e::PetscLogEvent, comm::MPI_Comm) 

# External Links
$(_doc_external("Sys/PetscLogEventSync"))
"""
function PetscLogEventSync(petsclib::PetscLibType, e::PetscLogEvent, comm::MPI_Comm) end

@for_petsc function PetscLogEventSync(petsclib::$UnionPetscLib, e::PetscLogEvent, comm::MPI_Comm )

    @chk ccall(
               (:PetscLogEventSync, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent, MPI_Comm),
               e, comm,
              )


	return nothing
end 

"""
	PetscLogStageGetPerfInfo(petsclib::PetscLibType,stage::PetscLogStage, info::PetscEventPerfInfo) 
Return the performance information about the given stage

No Fortran Support

Input Parameters:
- `stage` - The stage number or `PETSC_DETERMINE` for the current stage

Output Parameter:
- `info` - This structure is filled with the performance information

Level: intermediate

-seealso: [](ch_profiling), `PetscLogEventRegister()`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogStageGetPerfInfo"))
"""
function PetscLogStageGetPerfInfo(petsclib::PetscLibType, stage::PetscLogStage, info::PetscEventPerfInfo) end

@for_petsc function PetscLogStageGetPerfInfo(petsclib::$UnionPetscLib, stage::PetscLogStage, info::PetscEventPerfInfo )

    @chk ccall(
               (:PetscLogStageGetPerfInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, Ptr{PetscEventPerfInfo}),
               stage, info,
              )


	return nothing
end 

"""
	PetscLogEventGetPerfInfo(petsclib::PetscLibType,stage::PetscLogStage, event::PetscLogEvent, info::PetscEventPerfInfo) 
Return the performance information about the given event in the given stage

No Fortran Support

Input Parameters:
- `stage` - The stage number or `PETSC_DETERMINE` for the current stage
- `event` - The event number

Output Parameter:
- `info` - This structure is filled with the performance information

Level: intermediate

-seealso: [](ch_profiling), `PetscLogEventRegister()`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogEventGetPerfInfo"))
"""
function PetscLogEventGetPerfInfo(petsclib::PetscLibType, stage::PetscLogStage, event::PetscLogEvent, info::PetscEventPerfInfo) end

@for_petsc function PetscLogEventGetPerfInfo(petsclib::$UnionPetscLib, stage::PetscLogStage, event::PetscLogEvent, info::PetscEventPerfInfo )

    @chk ccall(
               (:PetscLogEventGetPerfInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogStage, PetscLogEvent, Ptr{PetscEventPerfInfo}),
               stage, event, info,
              )


	return nothing
end 

"""
	PetscLogEventSetDof(petsclib::PetscLibType,event::PetscLogEvent, n::PetscInt, dof::PetscLogDouble) 
Set the nth number of degrees of freedom of a numerical problem associated with this event

Not Collective

Input Parameters:
- `event` - The event id to log
- `n`     - The dof index, in [0, 8)
- `dof`   - The number of dofs

Options Database Key:
- `-log_view` - Activates log summary

Level: developer

-seealso: `PetscLogEventSetError()`, `PetscLogEventRegister()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogEventSetDof"))
"""
function PetscLogEventSetDof(petsclib::PetscLibType, event::PetscLogEvent, n::PetscInt, dof::PetscLogDouble) end

@for_petsc function PetscLogEventSetDof(petsclib::$UnionPetscLib, event::PetscLogEvent, n::$PetscInt, dof::PetscLogDouble )

    @chk ccall(
               (:PetscLogEventSetDof, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent, $PetscInt, PetscLogDouble),
               event, n, dof,
              )


	return nothing
end 

"""
	PetscLogEventSetError(petsclib::PetscLibType,event::PetscLogEvent, n::PetscInt, error::PetscLogDouble) 
Set the nth error associated with a numerical problem associated with this event

Not Collective

Input Parameters:
- `event` - The event id to log
- `n`     - The error index, in [0, 8)
- `error` - The error

Options Database Key:
- `-log_view` - Activates log summary

Level: developer

-seealso: `PetscLogEventSetDof()`, `PetscLogEventRegister()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogEventSetError"))
"""
function PetscLogEventSetError(petsclib::PetscLibType, event::PetscLogEvent, n::PetscInt, error::PetscLogDouble) end

@for_petsc function PetscLogEventSetError(petsclib::$UnionPetscLib, event::PetscLogEvent, n::$PetscInt, error::PetscLogDouble )

    @chk ccall(
               (:PetscLogEventSetError, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent, $PetscInt, PetscLogDouble),
               event, n, error,
              )


	return nothing
end 

"""
	PetscLogEventGetId(petsclib::PetscLibType,name::Vector{Cchar}, event::PetscLogEvent) 
Returns the event id when given the event name.

Not Collective

Input Parameter:
- `name` - The event name

Output Parameter:
- `event` - The event, or -1 if no event with that name exists

Level: intermediate

-seealso: [](ch_profiling), `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogStageGetId()`

# External Links
$(_doc_external("Sys/PetscLogEventGetId"))
"""
function PetscLogEventGetId(petsclib::PetscLibType, name::Vector{Cchar}, event::PetscLogEvent) end

@for_petsc function PetscLogEventGetId(petsclib::$UnionPetscLib, name::Vector{Cchar}, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogEventGetId, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscLogEvent}),
               name, event,
              )


	return nothing
end 

"""
	PetscLogEventGetName(petsclib::PetscLibType,event::PetscLogEvent, name::Vector{Cchar}) 
Returns the event name when given the event id.

Not Collective

Input Parameter:
- `event` - The event

Output Parameter:
- `name` - The event name

Level: intermediate

-seealso: [](ch_profiling), `PetscLogEventRegister()`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscPreLoadBegin()`, `PetscPreLoadEnd()`, `PetscPreLoadStage()`

# External Links
$(_doc_external("Sys/PetscLogEventGetName"))
"""
function PetscLogEventGetName(petsclib::PetscLibType, event::PetscLogEvent, name::Vector{Cchar}) end

@for_petsc function PetscLogEventGetName(petsclib::$UnionPetscLib, event::PetscLogEvent, name::Vector{Cchar} )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscLogEventGetName, $petsc_library),
               PetscErrorCode,
               (PetscLogEvent, Ptr{Ptr{Cchar}}),
               event, name_,
              )


	return nothing
end 

"""
	PetscLogEventsPause(petsclib::PetscLibType) 
Put event logging into "paused" mode: timers and counters for in

Not collective

Level: advanced

-seealso: [](ch_profiling), `PetscLogEventDeactivatePush()`, `PetscLogEventDeactivatePop()`, `PetscLogEventsResume()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogEventsPause"))
"""
function PetscLogEventsPause(petsclib::PetscLibType) end

@for_petsc function PetscLogEventsPause(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogEventsPause, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogEventsResume(petsclib::PetscLibType) 
Return logging to normal behavior after it was paused with `PetscLogEventsPause()`.

Not collective

Level: advanced

-seealso: [](ch_profiling), `PetscLogEventDeactivatePush()`, `PetscLogEventDeactivatePop()`, `PetscLogEventsPause()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogEventsResume"))
"""
function PetscLogEventsResume(petsclib::PetscLibType) end

@for_petsc function PetscLogEventsResume(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogEventsResume, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogObjectCreate(petsclib::PetscLibType,h::PetscObject) 

# External Links
$(_doc_external("Sys/PetscLogObjectCreate"))
"""
function PetscLogObjectCreate(petsclib::PetscLibType, h::PetscObject) end

@for_petsc function PetscLogObjectCreate(petsclib::$UnionPetscLib, h::PetscObject )

    @chk ccall(
               (:PetscLogObjectCreate, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               h,
              )


	return nothing
end 

"""
	PetscLogObjectDestroy(petsclib::PetscLibType,h::PetscObject) 

# External Links
$(_doc_external("Sys/PetscLogObjectDestroy"))
"""
function PetscLogObjectDestroy(petsclib::PetscLibType, h::PetscObject) end

@for_petsc function PetscLogObjectDestroy(petsclib::$UnionPetscLib, h::PetscObject )

    @chk ccall(
               (:PetscLogObjectDestroy, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               h,
              )


	return nothing
end 

"""
	PetscLogClassGetClassId(petsclib::PetscLibType,name::Vector{Cchar}, classid::PetscClassId) 
Returns the `PetscClassId` when given the class name.

Not Collective

Input Parameter:
- `name` - The class name

Output Parameter:
- `classid` - The `PetscClassId` id, or -1 if no class with that name exists

Level: intermediate

-seealso: [](ch_profiling), `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogStageGetId()`

# External Links
$(_doc_external("Sys/PetscLogClassGetClassId"))
"""
function PetscLogClassGetClassId(petsclib::PetscLibType, name::Vector{Cchar}, classid::PetscClassId) end

@for_petsc function PetscLogClassGetClassId(petsclib::$UnionPetscLib, name::Vector{Cchar}, classid::PetscClassId )

    @chk ccall(
               (:PetscLogClassGetClassId, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscClassId}),
               name, classid,
              )


	return nothing
end 

"""
	PetscLogClassIdGetName(petsclib::PetscLibType,classid::PetscClassId, name::Cchar) 
Returns a `PetscClassId`'s name.

Not Collective

Input Parameter:
- `classid` - A `PetscClassId`

Output Parameter:
- `name` - The class name

Level: intermediate

-seealso: [](ch_profiling), `PetscLogClassRegister()`, `PetscLogClassBegin()`, `PetscLogClassEnd()`, `PetscPreLoadBegin()`, `PetscPreLoadEnd()`, `PetscPreLoadClass()`

# External Links
$(_doc_external("Sys/PetscLogClassIdGetName"))
"""
function PetscLogClassIdGetName(petsclib::PetscLibType, classid::PetscClassId, name::Cchar) end

@for_petsc function PetscLogClassIdGetName(petsclib::$UnionPetscLib, classid::PetscClassId, name::Cchar )

    @chk ccall(
               (:PetscLogClassIdGetName, $petsc_library),
               PetscErrorCode,
               (PetscClassId, Cchar),
               classid, name,
              )


	return nothing
end 

"""
	PetscLogDump(petsclib::PetscLibType,sname::Vector{Cchar}) 
Dumps logs of objects to a file. This file is intended to
be read by bin/petscview. This program no longer exists.

Collective on `PETSC_COMM_WORLD`

Input Parameter:
- `sname` - an optional file name

-seealso: [](ch_profiling), `PetscLogDefaultBegin()`, `PetscLogView()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Sys/PetscLogDump"))
"""
function PetscLogDump(petsclib::PetscLibType, sname::Vector{Cchar}) end

@for_petsc function PetscLogDump(petsclib::$UnionPetscLib, sname::Vector{Cchar} )

    @chk ccall(
               (:PetscLogDump, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               sname,
              )


	return nothing
end 

"""
	PetscLogMPEDump(petsclib::PetscLibType,sname::Vector{Cchar}) 
Dumps the MPE logging info to file for later use with Jumpshot.

Collective on `PETSC_COMM_WORLD`

Input Parameter:
- `sname` - filename for the MPE logfile

Level: advanced

-seealso: [](ch_profiling), `PetscLogDump()`, `PetscLogMPEBegin()`

# External Links
$(_doc_external("Sys/PetscLogMPEDump"))
"""
function PetscLogMPEDump(petsclib::PetscLibType, sname::Vector{Cchar}) end

@for_petsc function PetscLogMPEDump(petsclib::$UnionPetscLib, sname::Vector{Cchar} )

    @chk ccall(
               (:PetscLogMPEDump, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               sname,
              )


	return nothing
end 

"""
	PetscLogView(petsclib::PetscLibType,viewer::PetscViewer) 
Prints a summary of the logging.

Collective

Input Parameter:
- `viewer` - an ASCII viewer

Options Database Keys:
- `-log_view [:filename]`                    - Prints summary of log information
- `-log_view :filename.py:ascii_info_detail` - Saves logging information from each process as a Python file
- `-log_view :filename.xml:ascii_xml`        - Saves a summary of the logging information in a nested format (see below for how to view it)
- `-log_view :filename.txt:ascii_flamegraph` - Saves logging information in a format suitable for visualising as a Flame Graph (see below for how to view it)
- `-log_view_memory`                         - Also display memory usage in each event
- `-log_view_gpu_time`                       - Also display time in each event for GPU kernels (Note this may slow the computation)
- `-log_all`                                 - Saves a file Log.rank for each MPI rank with details of each step of the computation
- `-log_trace [filename]`                    - Displays a trace of what each process is doing

Level: beginner

-seealso: [](ch_profiling), `PetscLogDefaultBegin()`, `PetscLogDump()`

# External Links
$(_doc_external("Sys/PetscLogView"))
"""
function PetscLogView(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscLogView(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscLogView, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscLogViewFromOptions(petsclib::PetscLibType) 
Processes command line options to determine if/how a `PetscLog` is to be viewed.

Collective on `PETSC_COMM_WORLD`

Level: developer

-seealso: [](ch_profiling), `PetscLogView()`

# External Links
$(_doc_external("Sys/PetscLogViewFromOptions"))
"""
function PetscLogViewFromOptions(petsclib::PetscLibType) end

@for_petsc function PetscLogViewFromOptions(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogViewFromOptions, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogSetThreshold(petsclib::PetscLibType,newThresh::PetscLogDouble, oldThresh::PetscLogDouble) 
Set the threshold time for logging the events; this is a percentage out of 100, so 1. means any event
that takes 1 or more percent of the time.

Logically Collective on `PETSC_COMM_WORLD`

Input Parameter:
- `newThresh` - the threshold to use

Output Parameter:
- `oldThresh` - the previously set threshold value

Options Database Keys:
- `-log_view :filename.xml:ascii_xml` - Prints an XML summary of flop and timing information to the file

-seealso: `PetscLogDump()`, `PetscLogView()`, `PetscLogTraceBegin()`, `PetscLogDefaultBegin()`,
`PetscLogNestedBegin()`

# External Links
$(_doc_external("Sys/PetscLogSetThreshold"))
"""
function PetscLogSetThreshold(petsclib::PetscLibType, newThresh::PetscLogDouble, oldThresh::PetscLogDouble) end

@for_petsc function PetscLogSetThreshold(petsclib::$UnionPetscLib, newThresh::PetscLogDouble, oldThresh::PetscLogDouble )

    @chk ccall(
               (:PetscLogSetThreshold, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble, Ptr{PetscLogDouble}),
               newThresh, oldThresh,
              )


	return nothing
end 

"""
	PetscGetFlops(petsclib::PetscLibType,flops::PetscLogDouble) 
Returns the number of flops used on this processor
since the program began.

Not Collective

Output Parameter:
- `flops` - number of floating point operations

Level: intermediate

-seealso: [](ch_profiling), `PetscLogGpuFlops()`, `PetscTime()`, `PetscLogFlops()`

# External Links
$(_doc_external("Sys/PetscGetFlops"))
"""
function PetscGetFlops(petsclib::PetscLibType, flops::PetscLogDouble) end

@for_petsc function PetscGetFlops(petsclib::$UnionPetscLib, flops::PetscLogDouble )

    @chk ccall(
               (:PetscGetFlops, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               flops,
              )


	return nothing
end 

"""
	PetscLogFlops(petsclib::PetscLibType,f::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogFlops"))
"""
function PetscLogFlops(petsclib::PetscLibType, f::PetscLogDouble) end

@for_petsc function PetscLogFlops(petsclib::$UnionPetscLib, f::PetscLogDouble )

    @chk ccall(
               (:PetscLogFlops, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               f,
              )


	return nothing
end 

"""
	PetscLogGpuTime(petsclib::PetscLibType) 
turn on the logging of GPU time for GPU kernels

Options Database Key:
- `-log_view_gpu_time` - provide the GPU times for all events in the `-log_view` output

Level: advanced

-seealso: [](ch_profiling), `PetscLogView()`, `PetscLogGpuFlops()`, `PetscLogGpuTimeEnd()`, `PetscLogGpuTimeBegin()`

# External Links
$(_doc_external("Sys/PetscLogGpuTime"))
"""
function PetscLogGpuTime(petsclib::PetscLibType) end

@for_petsc function PetscLogGpuTime(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogGpuTime, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogGpuTimeBegin(petsclib::PetscLibType) 
Start timer for device

Level: intermediate

-seealso: [](ch_profiling), `PetscLogView()`, `PetscLogGpuFlops()`, `PetscLogGpuTimeEnd()`, `PetscLogGpuTime()`

# External Links
$(_doc_external("Sys/PetscLogGpuTimeBegin"))
"""
function PetscLogGpuTimeBegin(petsclib::PetscLibType) end

@for_petsc function PetscLogGpuTimeBegin(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogGpuTimeBegin, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscLogGpuTimeEnd(petsclib::PetscLibType) 
Stop timer for device

Level: intermediate

-seealso: [](ch_profiling), `PetscLogView()`, `PetscLogGpuFlops()`, `PetscLogGpuTimeBegin()`

# External Links
$(_doc_external("Sys/PetscLogGpuTimeEnd"))
"""
function PetscLogGpuTimeEnd(petsclib::PetscLibType) end

@for_petsc function PetscLogGpuTimeEnd(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscLogGpuTimeEnd, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscClassIdRegister(petsclib::PetscLibType,name::Vector{Cchar}, oclass::PetscClassId) 
Registers a new class name for objects and logging operations in an application code.

Not Collective

Input Parameter:
- `name` - The class name

Output Parameter:
- `oclass` - The class id or classid

Level: developer

-seealso: [](ch_profiling), `PetscLogEventRegister()`

# External Links
$(_doc_external("Sys/PetscClassIdRegister"))
"""
function PetscClassIdRegister(petsclib::PetscLibType, name::Vector{Cchar}, oclass::PetscClassId) end

@for_petsc function PetscClassIdRegister(petsclib::$UnionPetscLib, name::Vector{Cchar}, oclass::PetscClassId )

    @chk ccall(
               (:PetscClassIdRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscClassId}),
               name, oclass,
              )


	return nothing
end 

"""
	keysout::Vector{PetscInt} = PetscParallelSortInt(petsclib::PetscLibType,mapin::PetscLayout, mapout::PetscLayout, keysin::Vector{PetscInt}) 
Globally sort a distributed array of integers

Collective

Input Parameters:
- `mapin`  - `PetscLayout` describing the distribution of the input keys
- `mapout` - `PetscLayout` describing the desired distribution of the output keys
- `keysin` - the pre-sorted array of integers

Output Parameter:
- `keysout` - the array in which the sorted integers will be stored.  If `mapin` == `mapout`, then `keysin` may be equal to `keysout`.

Level: developer

-seealso: `PetscSortInt()`, `PetscParallelSortedInt()`

# External Links
$(_doc_external("Vec/PetscParallelSortInt"))
"""
function PetscParallelSortInt(petsclib::PetscLibType, mapin::PetscLayout, mapout::PetscLayout, keysin::Vector{PetscInt}) end

@for_petsc function PetscParallelSortInt(petsclib::$UnionPetscLib, mapin::PetscLayout, mapout::PetscLayout, keysin::Vector{$PetscInt} )
	keysout = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscParallelSortInt, $petsc_library),
               PetscErrorCode,
               (PetscLayout, PetscLayout, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mapin, mapout, keysin, keysout,
              )


	return keysout
end 

"""
	PetscCommSplitReductionBegin(petsclib::PetscLibType,comm::MPI_Comm) 
Begin an asynchronous split

Collective but not synchronizing

Input Parameter:
- `comm` - communicator on which split reduction has been queued

Level: advanced

-seealso: `VecNormBegin()`, `VecNormEnd()`, `VecDotBegin()`, `VecDotEnd()`, `VecTDotBegin()`, `VecTDotEnd()`, `VecMDotBegin()`, `VecMDotEnd()`, `VecMTDotBegin()`, `VecMTDotEnd()`

# External Links
$(_doc_external("Vec/PetscCommSplitReductionBegin"))
"""
function PetscCommSplitReductionBegin(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscCommSplitReductionBegin(petsclib::$UnionPetscLib, comm::MPI_Comm )

    @chk ccall(
               (:PetscCommSplitReductionBegin, $petsc_library),
               PetscErrorCode,
               (MPI_Comm,),
               comm,
              )


	return nothing
end 

"""
	wv::PetscReal = PetscDTAltVApply(petsclib::PetscLibType,N::PetscInt, k::PetscInt, w::PetscReal, v::PetscReal) 
Apply an a k

Input Parameters:
- `N` - the dimension of the vector space, N >= 0
- `k` - the degree k of the k-form w, 0 <= k <= N
- `w` - a k-form, size [N choose k] (each degree of freedom of a k-form is associated with a subset of k coordinates of the N-dimensional vectors.
The degrees of freedom are ordered lexicographically by their associated subsets)
- `v` - a set of k vectors of size N, size [k x N], each vector stored contiguously

Output Parameter:
- `wv` 

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVApply"))
"""
function PetscDTAltVApply(petsclib::PetscLibType, N::PetscInt, k::PetscInt, w::PetscReal, v::PetscReal) end

@for_petsc function PetscDTAltVApply(petsclib::$UnionPetscLib, N::$PetscInt, k::$PetscInt, w::$PetscReal, v::$PetscReal )
	wv_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVApply, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, k, w, v, wv_,
              )

	wv = wv_[]

	return wv
end 

"""
	awedgeb::PetscReal = PetscDTAltVWedge(petsclib::PetscLibType,N::PetscInt, j::PetscInt, k::PetscInt, a::PetscReal, b::PetscReal) 
Compute the wedge product of a j

Input Parameters:
- `N` - the dimension of the vector space, N >= 0
- `j` - the degree j of the j-form a, 0 <= j <= N
- `k` - the degree k of the k-form b, 0 <= k <= N and 0 <= j+k <= N
- `a` - a j-form, size [N choose j]
- `b` - a k-form, size [N choose k]

Output Parameter:
- `awedgeb` - 
Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVWedgeMatrix()`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVWedge"))
"""
function PetscDTAltVWedge(petsclib::PetscLibType, N::PetscInt, j::PetscInt, k::PetscInt, a::PetscReal, b::PetscReal) end

@for_petsc function PetscDTAltVWedge(petsclib::$UnionPetscLib, N::$PetscInt, j::$PetscInt, k::$PetscInt, a::$PetscReal, b::$PetscReal )
	awedgeb_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVWedge, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, j, k, a, b, awedgeb_,
              )

	awedgeb = awedgeb_[]

	return awedgeb
end 

"""
	awedgeMat::PetscReal = PetscDTAltVWedgeMatrix(petsclib::PetscLibType,N::PetscInt, j::PetscInt, k::PetscInt, a::PetscReal) 
Compute the matrix defined by the wedge product with a given j

Input Parameters:
- `N` - the dimension of the vector space, N >= 0
- `j` - the degree j of the j-form a, 0 <= j <= N
- `k` - the degree k of the k-forms that (a wedge) will be applied to, 0 <= k <= N and 0 <= j+k <= N
- `a` - a j-form, size [N choose j]

Output Parameter:
- `awedgeMat` - (a wedge), an [(N choose j+k) x (N choose k)] matrix in row-major order, such that (a wedge) * b = a wedge b

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVWedgeMatrix"))
"""
function PetscDTAltVWedgeMatrix(petsclib::PetscLibType, N::PetscInt, j::PetscInt, k::PetscInt, a::PetscReal) end

@for_petsc function PetscDTAltVWedgeMatrix(petsclib::$UnionPetscLib, N::$PetscInt, j::$PetscInt, k::$PetscInt, a::$PetscReal )
	awedgeMat_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVWedgeMatrix, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, j, k, a, awedgeMat_,
              )

	awedgeMat = awedgeMat_[]

	return awedgeMat
end 

"""
	Lstarw::PetscReal = PetscDTAltVPullback(petsclib::PetscLibType,N::PetscInt, M::PetscInt, L::PetscReal, k::PetscInt, w::PetscReal) 
Compute the pullback of a k

Input Parameters:
- `N` - the dimension of the origin vector space of the linear transformation, M >= 0
- `M` - the dimension of the image vector space of the linear transformation, N >= 0
- `L` - a linear transformation, an [M x N] matrix in row-major format
- `k` - the *signed* degree k of the |k|-form w, -(min(M,N)) <= k <= min(M,N).  A negative form degree indicates that the pullback should be conjugated by
the Hodge star operator (see note).
- `w` - a |k|-form in the image space, size [M choose |k|]

Output Parameter:
- `Lstarw` - the pullback of w to a |k|-form in the origin space, size [N choose |k|]: (Lstarw)(v_1,...v_k) = w(L*v_1,...,L*v_k).

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVPullbackMatrix()`, `PetscDTAltVStar()`

# External Links
$(_doc_external("DM/PetscDTAltVPullback"))
"""
function PetscDTAltVPullback(petsclib::PetscLibType, N::PetscInt, M::PetscInt, L::PetscReal, k::PetscInt, w::PetscReal) end

@for_petsc function PetscDTAltVPullback(petsclib::$UnionPetscLib, N::$PetscInt, M::$PetscInt, L::$PetscReal, k::$PetscInt, w::$PetscReal )
	Lstarw_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVPullback, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, M, L, k, w, Lstarw_,
              )

	Lstarw = Lstarw_[]

	return Lstarw
end 

"""
	Lstar::PetscReal = PetscDTAltVPullbackMatrix(petsclib::PetscLibType,N::PetscInt, M::PetscInt, L::PetscReal, k::PetscInt) 
Compute the pullback matrix for k

Input Parameters:
- `N` - the dimension of the origin vector space of the linear transformation, N >= 0
- `M` - the dimension of the image vector space of the linear transformation, M >= 0
- `L` - a linear transformation, an [M x N] matrix in row-major format
- `k` - the *signed* degree k of the |k|-forms on which Lstar acts, -(min(M,N)) <= k <= min(M,N).
A negative form degree indicates that the pullback should be conjugated by the Hodge star operator (see note in `PetscDTAltvPullback()`)

Output Parameter:
- `Lstar` - the pullback matrix, an [(N choose |k|) x (M choose |k|)] matrix in row-major format such that Lstar * w = L^* w

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVPullback()`, `PetscDTAltVStar()`

# External Links
$(_doc_external("DM/PetscDTAltVPullbackMatrix"))
"""
function PetscDTAltVPullbackMatrix(petsclib::PetscLibType, N::PetscInt, M::PetscInt, L::PetscReal, k::PetscInt) end

@for_petsc function PetscDTAltVPullbackMatrix(petsclib::$UnionPetscLib, N::$PetscInt, M::$PetscInt, L::$PetscReal, k::$PetscInt )
	Lstar_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVPullbackMatrix, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}),
               N, M, L, k, Lstar_,
              )

	Lstar = Lstar_[]

	return Lstar
end 

"""
	wIntv::PetscReal = PetscDTAltVInterior(petsclib::PetscLibType,N::PetscInt, k::PetscInt, w::PetscReal, v::PetscReal) 
Compute the interior product of a k

Input Parameters:
- `N` - the dimension of the vector space, N >= 0
- `k` - the degree k of the k-form w, 0 <= k <= N
- `w` - a k-form, size [N choose k]
- `v` - an N dimensional vector

Output Parameter:
- `wIntv` - the (k-1)-form (w int v), size [N choose (k-1)]: (w int v) is defined by its action on (k-1) vectors {v_1, ..., v_{k-1}} as (w inv v)(v_1, ..., v_{k-1}) = w(v, v_1, ..., v_{k-1}).

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVInteriorMatrix()`, `PetscDTAltVInteriorPattern()`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVInterior"))
"""
function PetscDTAltVInterior(petsclib::PetscLibType, N::PetscInt, k::PetscInt, w::PetscReal, v::PetscReal) end

@for_petsc function PetscDTAltVInterior(petsclib::$UnionPetscLib, N::$PetscInt, k::$PetscInt, w::$PetscReal, v::$PetscReal )
	wIntv_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVInterior, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, k, w, v, wIntv_,
              )

	wIntv = wIntv_[]

	return wIntv
end 

"""
	intvMat::PetscReal = PetscDTAltVInteriorMatrix(petsclib::PetscLibType,N::PetscInt, k::PetscInt, v::PetscReal) 
Compute the matrix of the linear transformation induced on a k

Input Parameters:
- `N` - the dimension of the vector space, N >= 0
- `k` - the degree k of the k-forms on which intvMat acts, 0 <= k <= N
- `v` - an N dimensional vector

Output Parameter:
- `intvMat` - an [(N choose (k-1)) x (N choose k)] matrix, row-major: (intvMat) * w = (w int v)

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVInterior()`, `PetscDTAltVInteriorPattern()`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVInteriorMatrix"))
"""
function PetscDTAltVInteriorMatrix(petsclib::PetscLibType, N::PetscInt, k::PetscInt, v::PetscReal) end

@for_petsc function PetscDTAltVInteriorMatrix(petsclib::$UnionPetscLib, N::$PetscInt, k::$PetscInt, v::$PetscReal )
	intvMat_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVInteriorMatrix, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, k, v, intvMat_,
              )

	intvMat = intvMat_[]

	return intvMat
end 

"""
	PetscDTAltVInteriorPattern(petsclib::PetscLibType,N::PetscInt, k::PetscInt, ::PetscInt(indices) 
compute the sparsity and sign pattern of the interior product matrix computed in `PetscDTAltVInteriorMatrix()`

Input Parameters:
- `N` - the dimension of the vector space, N >= 0
- `k` - the degree of the k-forms on which `intvMat` from `PetscDTAltVInteriorMatrix()` acts,  0 le k le N .

Output Parameter:
- `indices` - The interior product matrix `intvMat` has dimensions [(N choose (k-1)) x (N choose k)] and has (N choose k) * k
non-zeros.  indices[i][0] and indices[i][1] are the row and column of a non-zero, and its value is equal to the vector
coordinate v[j] if indices[i][2] = j, or -v[j] if indices[i][2] = -(j+1)

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVInterior()`, `PetscDTAltVInteriorMatrix()`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVInteriorPattern"))
"""
function PetscDTAltVInteriorPattern(petsclib::PetscLibType, N::PetscInt, k::PetscInt, indices::Vector{PetscInt}) end

@for_petsc function PetscDTAltVInteriorPattern(petsclib::$UnionPetscLib, N::$PetscInt, k::$PetscInt, indices::Vector{PetscInt})

    @chk ccall(
               (:PetscDTAltVInteriorPattern, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}),
               N, k, indices,
              )


	return nothing
end 

"""
	starw::PetscReal = PetscDTAltVStar(petsclib::PetscLibType,N::PetscInt, k::PetscInt, pow::PetscInt, w::PetscReal) 
Apply a power of the Hodge star operator, which maps k

Input Parameters:
- `N`   - the dimension of the vector space, N >= 0
- `k`   - the degree k of the k-form w, 0 <= k <= N
- `pow` - the number of times to apply the Hodge star operator: pow < 0 indicates that the inverse of the Hodge star operator should be applied |pow| times.
- `w`   - a k-form, size [N choose k]

Output Parameter:
- `starw` - (star)^pow w

Level: intermediate

-seealso: `PetscDTAltV`, `PetscDTAltVPullback()`, `PetscDTAltVPullbackMatrix()`

# External Links
$(_doc_external("DM/PetscDTAltVStar"))
"""
function PetscDTAltVStar(petsclib::PetscLibType, N::PetscInt, k::PetscInt, pow::PetscInt, w::PetscReal) end

@for_petsc function PetscDTAltVStar(petsclib::$UnionPetscLib, N::$PetscInt, k::$PetscInt, pow::$PetscInt, w::$PetscReal )
	starw_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTAltVStar, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               N, k, pow, w, starw_,
              )

	starw = starw_[]

	return starw
end 

"""
	pdf::PetscProbFn,cdf::PetscProbFn,sampler::PetscProbFn = PetscProbCreateFromOptions(petsclib::PetscLibType,dim::PetscInt, prefix::Vector{Cchar}, name::Vector{Cchar}) 
Return the probability distribution specified by the arguments and options

Not Collective

Input Parameters:
- `dim`    - The dimension of sample points
- `prefix` - The options prefix, or `NULL`
- `name`   - The options database name for the probability distribution type

Output Parameters:
- `pdf`     - The PDF of this type, or `NULL`
- `cdf`     - The CDF of this type, or `NULL`
- `sampler` - The PDF sampler of this type, or `NULL`

Level: intermediate

-seealso: `PetscProbFn`, `PetscPDFMaxwellBoltzmann1D()`, `PetscPDFGaussian1D()`, `PetscPDFConstant1D()`

# External Links
$(_doc_external("DM/PetscProbCreateFromOptions"))
"""
function PetscProbCreateFromOptions(petsclib::PetscLibType, dim::PetscInt, prefix::Vector{Cchar}, name::Vector{Cchar}) end

@for_petsc function PetscProbCreateFromOptions(petsclib::$UnionPetscLib, dim::$PetscInt, prefix::Vector{Cchar}, name::Vector{Cchar} )
	pdf_ = Ref{PetscProbFn}()
	cdf_ = Ref{PetscProbFn}()
	sampler_ = Ref{PetscProbFn}()

    @chk ccall(
               (:PetscProbCreateFromOptions, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Cchar}, Ptr{Cchar}, PetscProbFn, PetscProbFn, PetscProbFn),
               dim, prefix, name, pdf_, cdf_, sampler_,
              )

	pdf = pdf_[]
	cdf = cdf_[]
	sampler = sampler_[]

	return pdf,cdf,sampler
end 

"""
	alpha::PetscReal = PetscProbComputeKSStatistic(petsclib::PetscLibType,v::PetscVec, cdf::PetscProbFn) 
Compute the Kolmogorov

Collective

Input Parameters:
- `v`   - The data vector, blocksize is the sample dimension
- `cdf` - The analytic CDF

Output Parameter:
- `alpha` - The KS statistic

Level: advanced

-seealso: `PetscProbComputeKSStatisticWeighted()`, `PetscProbComputeKSStatisticMagnitude()`, `PetscProbFn`

# External Links
$(_doc_external("DM/PetscProbComputeKSStatistic"))
"""
function PetscProbComputeKSStatistic(petsclib::PetscLibType, v::PetscVec, cdf::PetscProbFn) end

@for_petsc function PetscProbComputeKSStatistic(petsclib::$UnionPetscLib, v::PetscVec, cdf::PetscProbFn )
	alpha_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscProbComputeKSStatistic, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscProbFn}, Ptr{$PetscReal}),
               v, cdf, alpha_,
              )

	alpha = alpha_[]

	return alpha
end 

"""
	alpha::PetscReal = PetscProbComputeKSStatisticWeighted(petsclib::PetscLibType,v::PetscVec, w::PetscVec, cdf::PetscProbFn) 
Compute the Kolmogorov

Collective

Input Parameters:
- `v`   - The data vector, blocksize is the sample dimension
- `w`   - The vector of weights for each sample, instead of the default 1/n
- `cdf` - The analytic CDF

Output Parameter:
- `alpha` - The KS statistic

Level: advanced

-seealso: `PetscProbComputeKSStatistic()`, `PetscProbComputeKSStatisticMagnitude()`, `PetscProbFn`

# External Links
$(_doc_external("DM/PetscProbComputeKSStatisticWeighted"))
"""
function PetscProbComputeKSStatisticWeighted(petsclib::PetscLibType, v::PetscVec, w::PetscVec, cdf::PetscProbFn) end

@for_petsc function PetscProbComputeKSStatisticWeighted(petsclib::$UnionPetscLib, v::PetscVec, w::PetscVec, cdf::PetscProbFn )
	alpha_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscProbComputeKSStatisticWeighted, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{PetscProbFn}, Ptr{$PetscReal}),
               v, w, cdf, alpha_,
              )

	alpha = alpha_[]

	return alpha
end 

"""
	alpha::PetscReal = PetscProbComputeKSStatisticMagnitude(petsclib::PetscLibType,v::PetscVec, cdf::PetscProbFn) 
Compute the Kolmogorov

Collective

Input Parameters:
- `v`   - The data vector, blocksize is the sample dimension
- `cdf` - The analytic CDF

Output Parameter:
- `alpha` - The KS statistic

Level: advanced

-seealso: `PetscProbComputeKSStatistic()`, `PetscProbComputeKSStatisticWeighted()`, `PetscProbFn`

# External Links
$(_doc_external("DM/PetscProbComputeKSStatisticMagnitude"))
"""
function PetscProbComputeKSStatisticMagnitude(petsclib::PetscLibType, v::PetscVec, cdf::PetscProbFn) end

@for_petsc function PetscProbComputeKSStatisticMagnitude(petsclib::$UnionPetscLib, v::PetscVec, cdf::PetscProbFn )
	alpha_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscProbComputeKSStatisticMagnitude, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscProbFn}, Ptr{$PetscReal}),
               v, cdf, alpha_,
              )

	alpha = alpha_[]

	return alpha
end 

"""
	PetscFormKeySort(petsclib::PetscLibType,n::PetscInt, arr::Vector{PetscFormKey}) 
Sorts an array of `PetscFormKey` in place in increasing order.

Not Collective

Input Parameters:
- `n`   - number of values
- `arr` - array of `PetscFormKey`

Level: intermediate

-seealso: `PetscFormKey`, `PetscIntSortSemiOrdered()`, `PetscSortInt()`

# External Links
$(_doc_external("DM/PetscFormKeySort"))
"""
function PetscFormKeySort(petsclib::PetscLibType, n::PetscInt, arr::Vector{PetscFormKey}) end

@for_petsc function PetscFormKeySort(petsclib::$UnionPetscLib, n::$PetscInt, arr::Vector{PetscFormKey} )

    @chk ccall(
               (:PetscFormKeySort, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscFormKey}),
               n, arr,
              )


	return nothing
end 

"""
	norm::PetscReal = PetscDTJacobiNorm(petsclib::PetscLibType,alpha::PetscReal, beta::PetscReal, n::PetscInt) 
Compute the weighted L2 norm of a Jacobi polynomial.

Input Parameters:
- `alpha` - the left exponent > -1
- `beta`  - the right exponent > -1
- `n`     - the polynomial degree

Output Parameter:
- `norm` - the weighted L2 norm

Level: beginner

-seealso: `PetscQuadrature`, `PetscDTJacobiEval()`

# External Links
$(_doc_external("DM/PetscDTJacobiNorm"))
"""
function PetscDTJacobiNorm(petsclib::PetscLibType, alpha::PetscReal, beta::PetscReal, n::PetscInt) end

@for_petsc function PetscDTJacobiNorm(petsclib::$UnionPetscLib, alpha::$PetscReal, beta::$PetscReal, n::$PetscInt )
	norm_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTJacobiNorm, $petsc_library),
               PetscErrorCode,
               ($PetscReal, $PetscReal, $PetscInt, Ptr{$PetscReal}),
               alpha, beta, n, norm_,
              )

	norm = norm_[]

	return norm
end 

"""
	p::Vector{PetscReal} = PetscDTJacobiEvalJet(petsclib::PetscLibType,alpha::PetscReal, beta::PetscReal, npoints::PetscInt, points::Vector{PetscReal}, degree::PetscInt, k::PetscInt) 
Evaluate the jet (function and derivatives) of the Jacobi polynomials basis up to a given degree.

Input Parameters:
- `alpha`   - the left exponent of the weight
- `beta`    - the right exponetn of the weight
- `npoints` - the number of points to evaluate the polynomials at
- `points`  - [npoints] array of point coordinates
- `degree`  - the maximm degree polynomial space to evaluate, (degree + 1) will be evaluated total.
- `k`       - the maximum derivative to evaluate in the jet, (k + 1) will be evaluated total.

Output Parameter:
- `p` - an array containing the evaluations of the Jacobi polynomials's jets on the points.  the size is (degree + 1) x
(k + 1) x npoints, which also describes the order of the dimensions of this three-dimensional array: the first
(slowest varying) dimension is polynomial degree; the second dimension is derivative order; the third (fastest
varying) dimension is the index of the evaluation point.

Level: advanced

-seealso: `PetscDTJacobiEval()`, `PetscDTPKDEvalJet()`

# External Links
$(_doc_external("DM/PetscDTJacobiEvalJet"))
"""
function PetscDTJacobiEvalJet(petsclib::PetscLibType, alpha::PetscReal, beta::PetscReal, npoints::PetscInt, points::Vector{PetscReal}, degree::PetscInt, k::PetscInt) end

@for_petsc function PetscDTJacobiEvalJet(petsclib::$UnionPetscLib, alpha::$PetscReal, beta::$PetscReal, npoints::$PetscInt, points::Vector{$PetscReal}, degree::$PetscInt, k::$PetscInt )
	p = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTJacobiEvalJet, $petsc_library),
               PetscErrorCode,
               ($PetscReal, $PetscReal, $PetscInt, Ptr{$PetscReal}, $PetscInt, $PetscInt, Ptr{$PetscReal}),
               alpha, beta, npoints, points, degree, k, p,
              )


	return p
end 

"""
	B::Vector{PetscReal},D::Vector{PetscReal},D2::Vector{PetscReal} = PetscDTJacobiEval(petsclib::PetscLibType,npoints::PetscInt, alpha::PetscReal, beta::PetscReal, points::PetscReal, ndegree::PetscInt, degrees::PetscInt) 
evaluate Jacobi polynomials for the weight function (1.+x)^{\alpha} (1.
at points

Not Collective

Input Parameters:
- `npoints` - number of spatial points to evaluate at
- `alpha`   - the left exponent > -1
- `beta`    - the right exponent > -1
- `points`  - array of locations to evaluate at
- `ndegree` - number of basis degrees to evaluate
- `degrees` - sorted array of degrees to evaluate

Output Parameters:
- `B`  - row-oriented basis evaluation matrix B[point*ndegree + degree] (dimension npoints*ndegrees, allocated by caller) (or `NULL`)
- `D`  - row-oriented derivative evaluation matrix (or `NULL`)
- `D2` - row-oriented second derivative evaluation matrix (or `NULL`)

Level: intermediate

-seealso: `PetscDTGaussQuadrature()`, `PetscDTLegendreEval()`

# External Links
$(_doc_external("DM/PetscDTJacobiEval"))
"""
function PetscDTJacobiEval(petsclib::PetscLibType, npoints::PetscInt, alpha::PetscReal, beta::PetscReal, points::PetscReal, ndegree::PetscInt, degrees::PetscInt) end

@for_petsc function PetscDTJacobiEval(petsclib::$UnionPetscLib, npoints::$PetscInt, alpha::$PetscReal, beta::$PetscReal, points::$PetscReal, ndegree::$PetscInt, degrees::$PetscInt )
	B = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	D = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	D2 = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTJacobiEval, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscReal, $PetscReal, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               npoints, alpha, beta, points, ndegree, degrees, B, D, D2,
              )


	return B,D,D2
end 

"""
	B::Vector{PetscReal},D::Vector{PetscReal},D2::Vector{PetscReal} = PetscDTLegendreEval(petsclib::PetscLibType,npoints::PetscInt, points::PetscReal, ndegree::PetscInt, degrees::PetscInt) 
evaluate Legendre polynomials at points

Not Collective

Input Parameters:
- `npoints` - number of spatial points to evaluate at
- `points`  - array of locations to evaluate at
- `ndegree` - number of basis degrees to evaluate
- `degrees` - sorted array of degrees to evaluate

Output Parameters:
- `B`  - row-oriented basis evaluation matrix B[point*ndegree + degree] (dimension `npoints`*`ndegrees`, allocated by caller) (or `NULL`)
- `D`  - row-oriented derivative evaluation matrix (or `NULL`)
- `D2` - row-oriented second derivative evaluation matrix (or `NULL`)

Level: intermediate

-seealso: `PetscDTGaussQuadrature()`

# External Links
$(_doc_external("DM/PetscDTLegendreEval"))
"""
function PetscDTLegendreEval(petsclib::PetscLibType, npoints::PetscInt, points::PetscReal, ndegree::PetscInt, degrees::PetscInt) end

@for_petsc function PetscDTLegendreEval(petsclib::$UnionPetscLib, npoints::$PetscInt, points::$PetscReal, ndegree::$PetscInt, degrees::$PetscInt )
	B = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	D = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	D2 = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTLegendreEval, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               npoints, points, ndegree, degrees, B, D, D2,
              )


	return B,D,D2
end 

"""
	degtup::Vector{PetscInt} = PetscDTIndexToGradedOrder(petsclib::PetscLibType,len::PetscInt, index::PetscInt) 
convert an index into a tuple of monomial degrees in a graded order (that is, if the degree sum of tuple x is less than the degree sum of tuple y,
then the index of x is smaller than the index of y)

Input Parameters:
- `len`   - the desired length of the degree tuple
- `index` - the index to convert: should be >= 0

Output Parameter:
- `degtup` - filled with a tuple of degrees

Level: beginner

-seealso: `PetscDTGradedOrderToIndex()`

# External Links
$(_doc_external("DM/PetscDTIndexToGradedOrder"))
"""
function PetscDTIndexToGradedOrder(petsclib::PetscLibType, len::PetscInt, index::PetscInt) end

@for_petsc function PetscDTIndexToGradedOrder(petsclib::$UnionPetscLib, len::$PetscInt, index::$PetscInt )
	degtup = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTIndexToGradedOrder, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}),
               len, index, degtup,
              )


	return degtup
end 

"""
	index::PetscInt = PetscDTGradedOrderToIndex(petsclib::PetscLibType,len::PetscInt, degtup::Vector{PetscInt}) 
convert a tuple into an index in a graded order, the inverse of `PetscDTIndexToGradedOrder()`.

Input Parameters:
- `len`    - the length of the degree tuple
- `degtup` - tuple with this length

Output Parameter:
- `index` - index in graded order: >= 0

Level: beginner

-seealso: `PetscDTIndexToGradedOrder()`

# External Links
$(_doc_external("DM/PetscDTGradedOrderToIndex"))
"""
function PetscDTGradedOrderToIndex(petsclib::PetscLibType, len::PetscInt, degtup::Vector{PetscInt}) end

@for_petsc function PetscDTGradedOrderToIndex(petsclib::$UnionPetscLib, len::$PetscInt, degtup::Vector{$PetscInt} )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTGradedOrderToIndex, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               len, degtup, index_,
              )

	index = index_[]

	return index
end 

"""
	p::Vector{PetscReal} = PetscDTPKDEvalJet(petsclib::PetscLibType,dim::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, degree::PetscInt, k::PetscInt) 
Evaluate the jet (function and derivatives) of the Proriol
the space of polynomials up to a given degree.

Input Parameters:
- `dim`     - the number of variables in the multivariate polynomials
- `npoints` - the number of points to evaluate the polynomials at
- `points`  - [npoints x dim] array of point coordinates
- `degree`  - the degree (sum of degrees on the variables in a monomial) of the polynomial space to evaluate.  There are ((dim + degree) choose dim) polynomials in this space.
- `k`       - the maximum order partial derivative to evaluate in the jet.  There are (dim + k choose dim) partial derivatives
in the jet.  Choosing k = 0 means to evaluate just the function and no derivatives

Output Parameter:
- `p` - an array containing the evaluations of the PKD polynomials' jets on the points.  The size is ((dim + degree)
choose dim) x ((dim + k) choose dim) x npoints, which also describes the order of the dimensions of this
three-dimensional array: the first (slowest varying) dimension is basis function index; the second dimension is jet
index; the third (fastest varying) dimension is the index of the evaluation point.

Level: advanced

-seealso: `PetscDTGradedOrderToIndex()`, `PetscDTIndexToGradedOrder()`, `PetscDTJacobiEvalJet()`

# External Links
$(_doc_external("DM/PetscDTPKDEvalJet"))
"""
function PetscDTPKDEvalJet(petsclib::PetscLibType, dim::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, degree::PetscInt, k::PetscInt) end

@for_petsc function PetscDTPKDEvalJet(petsclib::$UnionPetscLib, dim::$PetscInt, npoints::$PetscInt, points::Vector{$PetscReal}, degree::$PetscInt, k::$PetscInt )
	p = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTPKDEvalJet, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, $PetscInt, Ptr{$PetscReal}),
               dim, npoints, points, degree, k, p,
              )


	return p
end 

"""
	size::PetscInt = PetscDTPTrimmedSize(petsclib::PetscLibType,dim::PetscInt, degree::PetscInt, formDegree::PetscInt) 
The size of the trimmed polynomial space of k
which can be evaluated in `PetscDTPTrimmedEvalJet()`.

Input Parameters:
- `dim`        - the number of variables in the multivariate polynomials
- `degree`     - the degree (sum of degrees on the variables in a monomial) of the trimmed polynomial space.
- `formDegree` - the degree of the form

Output Parameter:
- `size` - The number ((`dim` + `degree`) choose (`dim` + `formDegree`)) x ((`degree` + `formDegree` - 1) choose (`formDegree`))

Level: advanced

-seealso: `PetscDTPTrimmedEvalJet()`

# External Links
$(_doc_external("DM/PetscDTPTrimmedSize"))
"""
function PetscDTPTrimmedSize(petsclib::PetscLibType, dim::PetscInt, degree::PetscInt, formDegree::PetscInt) end

@for_petsc function PetscDTPTrimmedSize(petsclib::$UnionPetscLib, dim::$PetscInt, degree::$PetscInt, formDegree::$PetscInt )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTPTrimmedSize, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               dim, degree, formDegree, size_,
              )

	size = size_[]

	return size
end 

"""
	p::Vector{PetscReal} = PetscDTPTrimmedEvalJet(petsclib::PetscLibType,dim::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, degree::PetscInt, formDegree::PetscInt, jetDegree::PetscInt) 
Evaluate the jet (function and derivatives) of a basis of the trimmed polynomial k
a given degree.

Input Parameters:
- `dim`        - the number of variables in the multivariate polynomials
- `npoints`    - the number of points to evaluate the polynomials at
- `points`     - [npoints x dim] array of point coordinates
- `degree`     - the degree (sum of degrees on the variables in a monomial) of the trimmed polynomial space to evaluate.
There are ((dim + degree) choose (dim + formDegree)) x ((degree + formDegree - 1) choose (formDegree)) polynomials in this space.
(You can use `PetscDTPTrimmedSize()` to compute this size.)
- `formDegree` - the degree of the form
- `jetDegree`  - the maximum order partial derivative to evaluate in the jet.  There are ((dim + jetDegree) choose dim) partial derivatives
in the jet.  Choosing jetDegree = 0 means to evaluate just the function and no derivatives

Output Parameter:
- `p` - an array containing the evaluations of the PKD polynomials' jets on the points.

Level: advanced

-seealso: `PetscDTPKDEvalJet()`, `PetscDTPTrimmedSize()`

# External Links
$(_doc_external("DM/PetscDTPTrimmedEvalJet"))
"""
function PetscDTPTrimmedEvalJet(petsclib::PetscLibType, dim::PetscInt, npoints::PetscInt, points::Vector{PetscReal}, degree::PetscInt, formDegree::PetscInt, jetDegree::PetscInt) end

@for_petsc function PetscDTPTrimmedEvalJet(petsclib::$UnionPetscLib, dim::$PetscInt, npoints::$PetscInt, points::Vector{$PetscReal}, degree::$PetscInt, formDegree::$PetscInt, jetDegree::$PetscInt )
	p = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTPTrimmedEvalJet, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}),
               dim, npoints, points, degree, formDegree, jetDegree, p,
              )


	return p
end 

"""
	x::Vector{PetscReal},w::Vector{PetscReal} = PetscDTGaussJacobiQuadrature(petsclib::PetscLibType,npoints::PetscInt, a::PetscReal, b::PetscReal, alpha::PetscReal, beta::PetscReal) 
quadrature for the interval [a, b] with the weight function
(x-a)^\alpha (x-b)^\beta.

Not Collective

Input Parameters:
- `npoints` - the number of points in the quadrature rule
- `a`       - the left endpoint of the interval
- `b`       - the right endpoint of the interval
- `alpha`   - the left exponent
- `beta`    - the right exponent

Output Parameters:
- `x` - array of length `npoints`, the locations of the quadrature points
- `w` - array of length `npoints`, the weights of the quadrature points

Level: intermediate

-seealso: `PetscDTGaussQuadrature()`

# External Links
$(_doc_external("DM/PetscDTGaussJacobiQuadrature"))
"""
function PetscDTGaussJacobiQuadrature(petsclib::PetscLibType, npoints::PetscInt, a::PetscReal, b::PetscReal, alpha::PetscReal, beta::PetscReal) end

@for_petsc function PetscDTGaussJacobiQuadrature(petsclib::$UnionPetscLib, npoints::$PetscInt, a::$PetscReal, b::$PetscReal, alpha::$PetscReal, beta::$PetscReal )
	x = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	w = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTGaussJacobiQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
               npoints, a, b, alpha, beta, x, w,
              )


	return x,w
end 

"""
	x::Vector{PetscReal},w::Vector{PetscReal} = PetscDTGaussLobattoJacobiQuadrature(petsclib::PetscLibType,npoints::PetscInt, a::PetscReal, b::PetscReal, alpha::PetscReal, beta::PetscReal) 
quadrature for the interval [a, b] with the weight function
(x-a)^\alpha (x-b)^\beta, with endpoints `a` and `b` included as quadrature points.

Not Collective

Input Parameters:
- `npoints` - the number of points in the quadrature rule
- `a`       - the left endpoint of the interval
- `b`       - the right endpoint of the interval
- `alpha`   - the left exponent
- `beta`    - the right exponent

Output Parameters:
- `x` - array of length `npoints`, the locations of the quadrature points
- `w` - array of length `npoints`, the weights of the quadrature points

Level: intermediate

-seealso: `PetscDTGaussJacobiQuadrature()`

# External Links
$(_doc_external("DM/PetscDTGaussLobattoJacobiQuadrature"))
"""
function PetscDTGaussLobattoJacobiQuadrature(petsclib::PetscLibType, npoints::PetscInt, a::PetscReal, b::PetscReal, alpha::PetscReal, beta::PetscReal) end

@for_petsc function PetscDTGaussLobattoJacobiQuadrature(petsclib::$UnionPetscLib, npoints::$PetscInt, a::$PetscReal, b::$PetscReal, alpha::$PetscReal, beta::$PetscReal )
	x = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	w = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTGaussLobattoJacobiQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
               npoints, a, b, alpha, beta, x, w,
              )


	return x,w
end 

"""
	x::Vector{PetscReal},w::Vector{PetscReal} = PetscDTGaussQuadrature(petsclib::PetscLibType,npoints::PetscInt, a::PetscReal, b::PetscReal) 
create Gauss

Not Collective

Input Parameters:
- `npoints` - number of points
- `a`       - left end of interval (often-1)
- `b`       - right end of interval (often +1)

Output Parameters:
- `x` - quadrature points
- `w` - quadrature weights

Level: intermediate

-seealso: `PetscDTLegendreEval()`, `PetscDTGaussJacobiQuadrature()`

# External Links
$(_doc_external("DM/PetscDTGaussQuadrature"))
"""
function PetscDTGaussQuadrature(petsclib::PetscLibType, npoints::PetscInt, a::PetscReal, b::PetscReal) end

@for_petsc function PetscDTGaussQuadrature(petsclib::$UnionPetscLib, npoints::$PetscInt, a::$PetscReal, b::$PetscReal )
	x = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	w = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTGaussQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscReal, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
               npoints, a, b, x, w,
              )


	return x,w
end 

"""
	x::Vector{PetscReal},w::Vector{PetscReal} = PetscDTGaussLobattoLegendreQuadrature(petsclib::PetscLibType,npoints::PetscInt, type::PetscGaussLobattoLegendreCreateType) 
creates a set of the locations and weights of the Gauss
nodes of a given size on the domain [-1,1]

Not Collective

Input Parameters:
- `npoints` - number of grid nodes
- `type`    - `PETSCGAUSSLOBATTOLEGENDRE_VIA_LINEAR_ALGEBRA` or `PETSCGAUSSLOBATTOLEGENDRE_VIA_NEWTON`

Output Parameters:
- `x` - quadrature points, pass in an array of length `npoints`
- `w` - quadrature weights, pass in an array of length `npoints`

Level: intermediate

-seealso: `PetscDTGaussQuadrature()`, `PetscGaussLobattoLegendreCreateType`


# External Links
$(_doc_external("DM/PetscDTGaussLobattoLegendreQuadrature"))
"""
function PetscDTGaussLobattoLegendreQuadrature(petsclib::PetscLibType, npoints::PetscInt, type::PetscGaussLobattoLegendreCreateType) end

@for_petsc function PetscDTGaussLobattoLegendreQuadrature(petsclib::$UnionPetscLib, npoints::$PetscInt, type::PetscGaussLobattoLegendreCreateType )
	x = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	w = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTGaussLobattoLegendreQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, PetscGaussLobattoLegendreCreateType, Ptr{$PetscReal}, Ptr{$PetscReal}),
               npoints, type, x, w,
              )


	return x,w
end 

"""
	PetscDTGaussTensorQuadrature(petsclib::PetscLibType,dim::PetscInt, Nc::PetscInt, npoints::PetscInt, a::PetscReal, b::PetscReal, q::PetscQuadrature) 
creates a tensor

Not Collective

Input Parameters:
- `dim`     - The spatial dimension
- `Nc`      - The number of components
- `npoints` - number of points in one dimension
- `a`       - left end of interval (often-1)
- `b`       - right end of interval (often +1)

Output Parameter:
- `q` - A `PetscQuadrature` object

Level: intermediate

-seealso: `PetscDTGaussQuadrature()`, `PetscDTLegendreEval()`

# External Links
$(_doc_external("DM/PetscDTGaussTensorQuadrature"))
"""
function PetscDTGaussTensorQuadrature(petsclib::PetscLibType, dim::PetscInt, Nc::PetscInt, npoints::PetscInt, a::PetscReal, b::PetscReal, q::PetscQuadrature) end

@for_petsc function PetscDTGaussTensorQuadrature(petsclib::$UnionPetscLib, dim::$PetscInt, Nc::$PetscInt, npoints::$PetscInt, a::$PetscReal, b::$PetscReal, q::PetscQuadrature )

    @chk ccall(
               (:PetscDTGaussTensorQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, $PetscReal, $PetscReal, Ptr{PetscQuadrature}),
               dim, Nc, npoints, a, b, q,
              )


	return nothing
end 

"""
	PetscDTStroudConicalQuadrature(petsclib::PetscLibType,dim::PetscInt, Nc::PetscInt, npoints::PetscInt, a::PetscReal, b::PetscReal, q::PetscQuadrature) 
create Stroud conical quadrature for a simplex {cite}`karniadakis2005spectral`

Not Collective

Input Parameters:
- `dim`     - The simplex dimension
- `Nc`      - The number of components
- `npoints` - The number of points in one dimension
- `a`       - left end of interval (often-1)
- `b`       - right end of interval (often +1)

Output Parameter:
- `q` - A `PetscQuadrature` object

Level: intermediate

-seealso: `PetscDTGaussTensorQuadrature()`, `PetscDTGaussQuadrature()`

# External Links
$(_doc_external("DM/PetscDTStroudConicalQuadrature"))
"""
function PetscDTStroudConicalQuadrature(petsclib::PetscLibType, dim::PetscInt, Nc::PetscInt, npoints::PetscInt, a::PetscReal, b::PetscReal, q::PetscQuadrature) end

@for_petsc function PetscDTStroudConicalQuadrature(petsclib::$UnionPetscLib, dim::$PetscInt, Nc::$PetscInt, npoints::$PetscInt, a::$PetscReal, b::$PetscReal, q::PetscQuadrature )

    @chk ccall(
               (:PetscDTStroudConicalQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, $PetscReal, $PetscReal, Ptr{PetscQuadrature}),
               dim, Nc, npoints, a, b, q,
              )


	return nothing
end 

"""
	PetscDTSimplexQuadrature(petsclib::PetscLibType,dim::PetscInt, degree::PetscInt, type::PetscDTSimplexQuadratureType, quad::PetscQuadrature) 
Create a quadrature rule for a simplex that exactly integrates polynomials up to a given degree.

Not Collective

Input Parameters:
- `dim`    - The spatial dimension of the simplex (1 = segment, 2 = triangle, 3 = tetrahedron)
- `degree` - The largest polynomial degree that is required to be integrated exactly
- `type`   - `PetscDTSimplexQuadratureType` indicating the type of quadrature rule

Output Parameter:
- `quad` - A `PetscQuadrature` object for integration over the biunit simplex

Level: intermediate

-seealso: `PetscDTSimplexQuadratureType`, `PetscDTGaussQuadrature()`, `PetscDTStroudCononicalQuadrature()`, `PetscQuadrature`

# External Links
$(_doc_external("DM/PetscDTSimplexQuadrature"))
"""
function PetscDTSimplexQuadrature(petsclib::PetscLibType, dim::PetscInt, degree::PetscInt, type::PetscDTSimplexQuadratureType, quad::PetscQuadrature) end

@for_petsc function PetscDTSimplexQuadrature(petsclib::$UnionPetscLib, dim::$PetscInt, degree::$PetscInt, type::PetscDTSimplexQuadratureType, quad::PetscQuadrature )

    @chk ccall(
               (:PetscDTSimplexQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, PetscDTSimplexQuadratureType, Ptr{PetscQuadrature}),
               dim, degree, type, quad,
              )


	return nothing
end 

"""
	PetscDTTanhSinhTensorQuadrature(petsclib::PetscLibType,dim::PetscInt, level::PetscInt, a::PetscReal, b::PetscReal, q::PetscQuadrature) 
create tanh

Not Collective

Input Parameters:
- `dim`   - The cell dimension
- `level` - The number of points in one dimension, 2^l
- `a`     - left end of interval (often-1)
- `b`     - right end of interval (often +1)

Output Parameter:
- `q` - A `PetscQuadrature` object

Level: intermediate

-seealso: `PetscDTGaussTensorQuadrature()`, `PetscQuadrature`

# External Links
$(_doc_external("DM/PetscDTTanhSinhTensorQuadrature"))
"""
function PetscDTTanhSinhTensorQuadrature(petsclib::PetscLibType, dim::PetscInt, level::PetscInt, a::PetscReal, b::PetscReal, q::PetscQuadrature) end

@for_petsc function PetscDTTanhSinhTensorQuadrature(petsclib::$UnionPetscLib, dim::$PetscInt, level::$PetscInt, a::$PetscReal, b::$PetscReal, q::PetscQuadrature )

    @chk ccall(
               (:PetscDTTanhSinhTensorQuadrature, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscReal, $PetscReal, Ptr{PetscQuadrature}),
               dim, level, a, b, q,
              )


	return nothing
end 

"""
	sol::PetscReal = PetscDTTanhSinhIntegrate(petsclib::PetscLibType,func::external, a::PetscReal, b::PetscReal, digits::PetscInt, ctx::Cvoid) 

# External Links
$(_doc_external("DM/PetscDTTanhSinhIntegrate"))
"""
function PetscDTTanhSinhIntegrate(petsclib::PetscLibType, func::external, a::PetscReal, b::PetscReal, digits::PetscInt, ctx::Cvoid) end

@for_petsc function PetscDTTanhSinhIntegrate(petsclib::$UnionPetscLib, func::external, a::$PetscReal, b::$PetscReal, digits::$PetscInt, ctx::Cvoid )
	sol_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTTanhSinhIntegrate, $petsc_library),
               PetscErrorCode,
               (external, $PetscReal, $PetscReal, $PetscInt, Ptr{Cvoid}, Ptr{$PetscReal}),
               func, a, b, digits, ctx, sol_,
              )

	sol = sol_[]

	return sol
end 

"""
	sol::PetscReal = PetscDTTanhSinhIntegrateMPFR(petsclib::PetscLibType,func::external, a::PetscReal, b::PetscReal, digits::PetscInt, ctx::Cvoid) 

# External Links
$(_doc_external("DM/PetscDTTanhSinhIntegrateMPFR"))
"""
function PetscDTTanhSinhIntegrateMPFR(petsclib::PetscLibType, func::external, a::PetscReal, b::PetscReal, digits::PetscInt, ctx::Cvoid) end

@for_petsc function PetscDTTanhSinhIntegrateMPFR(petsclib::$UnionPetscLib, func::external, a::$PetscReal, b::$PetscReal, digits::$PetscInt, ctx::Cvoid )
	sol_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTTanhSinhIntegrateMPFR, $petsc_library),
               PetscErrorCode,
               (external, $PetscReal, $PetscReal, $PetscInt, Ptr{Cvoid}, Ptr{$PetscReal}),
               func, a, b, digits, ctx, sol_,
              )

	sol = sol_[]

	return sol
end 

"""
	q::PetscQuadrature = PetscDTTensorQuadratureCreate(petsclib::PetscLibType,q1::PetscQuadrature, q2::PetscQuadrature) 
create the tensor product quadrature from two lower

Not Collective

Input Parameters:
- `q1` - The first quadrature
- `q2` - The second quadrature

Output Parameter:
- `q` - A `PetscQuadrature` object

Level: intermediate

-seealso: `PetscQuadrature`, `PetscDTGaussTensorQuadrature()`

# External Links
$(_doc_external("DM/PetscDTTensorQuadratureCreate"))
"""
function PetscDTTensorQuadratureCreate(petsclib::PetscLibType, q1::PetscQuadrature, q2::PetscQuadrature) end

@for_petsc function PetscDTTensorQuadratureCreate(petsclib::$UnionPetscLib, q1::PetscQuadrature, q2::PetscQuadrature )
	q_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:PetscDTTensorQuadratureCreate, $petsc_library),
               PetscErrorCode,
               (PetscQuadrature, PetscQuadrature, Ptr{PetscQuadrature}),
               q1, q2, q_,
              )

	q = q_[]

	return q
end 

"""
	R::Vector{PetscReal} = PetscDTReconstructPoly(petsclib::PetscLibType,degree::PetscInt, nsource::PetscInt, sourcex::Vector{PetscReal}, ntarget::PetscInt, targetx::Vector{PetscReal}) 
create matrix representing polynomial reconstruction using cell intervals and evaluation at target intervals

Not Collective

Input Parameters:
- `degree`  - degree of reconstruction polynomial
- `nsource` - number of source intervals
- `sourcex` - sorted coordinates of source cell boundaries (length `nsource`+1)
- `ntarget` - number of target intervals
- `targetx` - sorted coordinates of target cell boundaries (length `ntarget`+1)

Output Parameter:
- `R` - reconstruction matrix, utarget = sum_s R[t*nsource+s] * usource[s]

Level: advanced

-seealso: `PetscDTLegendreEval()`

# External Links
$(_doc_external("DM/PetscDTReconstructPoly"))
"""
function PetscDTReconstructPoly(petsclib::PetscLibType, degree::PetscInt, nsource::PetscInt, sourcex::Vector{PetscReal}, ntarget::PetscInt, targetx::Vector{PetscReal}) end

@for_petsc function PetscDTReconstructPoly(petsclib::$UnionPetscLib, degree::$PetscInt, nsource::$PetscInt, sourcex::Vector{$PetscReal}, ntarget::$PetscInt, targetx::Vector{$PetscReal} )
	R = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTReconstructPoly, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               degree, nsource, sourcex, ntarget, targetx, R,
              )


	return R
end 

"""
	in::PetscReal = PetscGaussLobattoLegendreIntegrate(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, f::Vector{PetscReal}) 
Compute the L2 integral of a function on the GLL points

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes
- `weights` - the GLL weights
- `f`       - the function values at the nodes

Output Parameter:
- `in` - the value of the integral

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreIntegrate"))
"""
function PetscGaussLobattoLegendreIntegrate(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, f::Vector{PetscReal}) end

@for_petsc function PetscGaussLobattoLegendreIntegrate(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal}, f::Vector{$PetscReal} )
	in_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscGaussLobattoLegendreIntegrate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               n, nodes, weights, f, in_,
              )

	in = in_[]

	return in
end 

"""
	AA::PetscReal = PetscGaussLobattoLegendreElementLaplacianCreate(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}) 
computes the Laplacian for a single 1d GLL element

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes, of length `n`
- `weights` - the GLL weights, of length `n`

Output Parameter:
- `AA` - the stiffness element, of size `n` by `n`

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`, `PetscGaussLobattoLegendreElementLaplacianDestroy()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementLaplacianCreate"))
"""
function PetscGaussLobattoLegendreElementLaplacianCreate(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}) end

@for_petsc function PetscGaussLobattoLegendreElementLaplacianCreate(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal} )
	AA_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscGaussLobattoLegendreElementLaplacianCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal),
               n, nodes, weights, AA_,
              )

	AA = AA_[]

	return AA
end 

"""
	PetscGaussLobattoLegendreElementLaplacianDestroy(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, AA::PetscReal) 
frees the Laplacian for a single 1d GLL element created with `PetscGaussLobattoLegendreElementLaplacianCreate()`

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes, ignored
- `weights` - the GLL weightss, ignored
- `AA`      - the stiffness element from `PetscGaussLobattoLegendreElementLaplacianCreate()`

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`, `PetscGaussLobattoLegendreElementLaplacianCreate()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementLaplacianDestroy"))
"""
function PetscGaussLobattoLegendreElementLaplacianDestroy(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, AA::PetscReal) end

@for_petsc function PetscGaussLobattoLegendreElementLaplacianDestroy(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal}, AA::$PetscReal )

    @chk ccall(
               (:PetscGaussLobattoLegendreElementLaplacianDestroy, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal),
               n, nodes, weights, AA,
              )


	return nothing
end 

"""
	AA::PetscReal,AAT::PetscReal = PetscGaussLobattoLegendreElementGradientCreate(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}) 
computes the gradient for a single 1d GLL element

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes, of length `n`
- `weights` - the GLL weights, of length `n`

Output Parameters:
- `AA`  - the stiffness element, of dimension `n` by `n`
- `AAT` - the transpose of AA (pass in `NULL` if you do not need this array), of dimension `n` by `n`

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`, `PetscGaussLobattoLegendreElementLaplacianDestroy()`, `PetscGaussLobattoLegendreElementGradientDestroy()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementGradientCreate"))
"""
function PetscGaussLobattoLegendreElementGradientCreate(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}) end

@for_petsc function PetscGaussLobattoLegendreElementGradientCreate(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal} )
	AA_ = Ref{$PetscReal}()
	AAT_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscGaussLobattoLegendreElementGradientCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal, $PetscReal),
               n, nodes, weights, AA_, AAT_,
              )

	AA = AA_[]
	AAT = AAT_[]

	return AA,AAT
end 

"""
	PetscGaussLobattoLegendreElementGradientDestroy(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, AA::PetscReal, AAT::PetscReal) 
frees the gradient for a single 1d GLL element obtained with `PetscGaussLobattoLegendreElementGradientCreate()`

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes, ignored
- `weights` - the GLL weights, ignored
- `AA`      - the stiffness element obtained with `PetscGaussLobattoLegendreElementGradientCreate()`
- `AAT`     - the transpose of the element obtained with `PetscGaussLobattoLegendreElementGradientCreate()`

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`, `PetscGaussLobattoLegendreElementLaplacianCreate()`, `PetscGaussLobattoLegendreElementAdvectionCreate()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementGradientDestroy"))
"""
function PetscGaussLobattoLegendreElementGradientDestroy(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, AA::PetscReal, AAT::PetscReal) end

@for_petsc function PetscGaussLobattoLegendreElementGradientDestroy(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal}, AA::$PetscReal, AAT::$PetscReal )

    @chk ccall(
               (:PetscGaussLobattoLegendreElementGradientDestroy, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal, $PetscReal),
               n, nodes, weights, AA, AAT,
              )


	return nothing
end 

"""
	AA::PetscReal = PetscGaussLobattoLegendreElementAdvectionCreate(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}) 
computes the advection operator for a single 1d GLL element

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes, of length `n`
- `weights` - the GLL weights, of length `n`

Output Parameter:
- `AA` - the stiffness element, of dimension `n` by `n`

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`, `PetscGaussLobattoLegendreElementLaplacianCreate()`, `PetscGaussLobattoLegendreElementAdvectionDestroy()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementAdvectionCreate"))
"""
function PetscGaussLobattoLegendreElementAdvectionCreate(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}) end

@for_petsc function PetscGaussLobattoLegendreElementAdvectionCreate(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal} )
	AA_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscGaussLobattoLegendreElementAdvectionCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal),
               n, nodes, weights, AA_,
              )

	AA = AA_[]

	return AA
end 

"""
	PetscGaussLobattoLegendreElementAdvectionDestroy(petsclib::PetscLibType,n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, AA::PetscReal) 
frees the advection stiffness for a single 1d GLL element created with `PetscGaussLobattoLegendreElementAdvectionCreate()`

Not Collective

Input Parameters:
- `n`       - the number of GLL nodes
- `nodes`   - the GLL nodes, ignored
- `weights` - the GLL weights, ignored
- `AA`      - advection obtained with `PetscGaussLobattoLegendreElementAdvectionCreate()`

Level: beginner

-seealso: `PetscDTGaussLobattoLegendreQuadrature()`, `PetscGaussLobattoLegendreElementAdvectionCreate()`

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementAdvectionDestroy"))
"""
function PetscGaussLobattoLegendreElementAdvectionDestroy(petsclib::PetscLibType, n::PetscInt, nodes::Vector{PetscReal}, weights::Vector{PetscReal}, AA::PetscReal) end

@for_petsc function PetscGaussLobattoLegendreElementAdvectionDestroy(petsclib::$UnionPetscLib, n::$PetscInt, nodes::Vector{$PetscReal}, weights::Vector{$PetscReal}, AA::$PetscReal )

    @chk ccall(
               (:PetscGaussLobattoLegendreElementAdvectionDestroy, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal),
               n, nodes, weights, AA,
              )


	return nothing
end 

"""
	nodes::PetscReal,weights::PetscReal = PetscGaussLobattoLegendreElementMassCreate(petsclib::PetscLibType,n::PetscInt, AA::PetscReal) 

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementMassCreate"))
"""
function PetscGaussLobattoLegendreElementMassCreate(petsclib::PetscLibType, n::PetscInt, AA::PetscReal) end

@for_petsc function PetscGaussLobattoLegendreElementMassCreate(petsclib::$UnionPetscLib, n::$PetscInt, AA::$PetscReal )
	nodes_ = Ref{$PetscReal}()
	weights_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscGaussLobattoLegendreElementMassCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal),
               n, nodes_, weights_, AA,
              )

	nodes = nodes_[]
	weights = weights_[]

	return nodes,weights
end 

"""
	nodes::PetscReal,weights::PetscReal = PetscGaussLobattoLegendreElementMassDestroy(petsclib::PetscLibType,n::PetscInt, AA::PetscReal) 

# External Links
$(_doc_external("DM/PetscGaussLobattoLegendreElementMassDestroy"))
"""
function PetscGaussLobattoLegendreElementMassDestroy(petsclib::PetscLibType, n::PetscInt, AA::PetscReal) end

@for_petsc function PetscGaussLobattoLegendreElementMassDestroy(petsclib::$UnionPetscLib, n::$PetscInt, AA::$PetscReal )
	nodes_ = Ref{$PetscReal}()
	weights_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscGaussLobattoLegendreElementMassDestroy, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal),
               n, nodes_, weights_, AA,
              )

	nodes = nodes_[]
	weights = weights_[]

	return nodes,weights
end 

"""
	coord::Vector{PetscInt} = PetscDTIndexToBary(petsclib::PetscLibType,len::PetscInt, sum::PetscInt, index::PetscInt) 
convert an index into a barycentric coordinate.

Input Parameters:
- `len`   - the desired length of the barycentric tuple (usually 1 more than the dimension it represents, so a barycentric coordinate in a triangle has length 3)
- `sum`   - the value that the sum of the barycentric coordinates (which will be non-negative integers) should sum to
- `index` - the index to convert: should be >= 0 and < Binomial(len - 1 + sum, sum)

Output Parameter:
- `coord` - will be filled with the barycentric coordinate, of length `n`

Level: beginner

-seealso: `PetscDTBaryToIndex()`

# External Links
$(_doc_external("DM/PetscDTIndexToBary"))
"""
function PetscDTIndexToBary(petsclib::PetscLibType, len::PetscInt, sum::PetscInt, index::PetscInt) end

@for_petsc function PetscDTIndexToBary(petsclib::$UnionPetscLib, len::$PetscInt, sum::$PetscInt, index::$PetscInt )
	coord = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscDTIndexToBary, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               len, sum, index, coord,
              )


	return coord
end 

"""
	index::PetscInt = PetscDTBaryToIndex(petsclib::PetscLibType,len::PetscInt, sum::PetscInt, coord::Vector{PetscInt}) 
convert a barycentric coordinate to an index

Input Parameters:
- `len`   - the desired length of the barycentric tuple (usually 1 more than the dimension it represents, so a barycentric coordinate in a triangle has length 3)
- `sum`   - the value that the sum of the barycentric coordinates (which will be non-negative integers) should sum to
- `coord` - a barycentric coordinate with the given length `len` and `sum`

Output Parameter:
- `index` - the unique index for the coordinate, >= 0 and < Binomial(len - 1 + sum, sum)

Level: beginner

-seealso: `PetscDTIndexToBary`

# External Links
$(_doc_external("DM/PetscDTBaryToIndex"))
"""
function PetscDTBaryToIndex(petsclib::PetscLibType, len::PetscInt, sum::PetscInt, coord::Vector{PetscInt}) end

@for_petsc function PetscDTBaryToIndex(petsclib::$UnionPetscLib, len::$PetscInt, sum::$PetscInt, coord::Vector{$PetscInt} )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTBaryToIndex, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               len, sum, coord, index_,
              )

	index = index_[]

	return index
end 

"""
	q::PetscQuadrature,fq::PetscQuadrature = PetscDTCreateQuadratureByCell(petsclib::PetscLibType,ct::DMPolytopeType, qorder::PetscInt, qtype::PetscDTSimplexQuadratureType) 
Create default quadrature for a given cell

Not collective

Input Parameters:
- `ct`     - The integration domain
- `qorder` - The desired quadrature order
- `qtype`  - The type of simplex quadrature, or PETSCDTSIMPLEXQUAD_DEFAULT

Output Parameters:
- `q`  - The cell quadrature
- `fq` - The face quadrature

Level: developer

-seealso: `PetscDTCreateDefaultQuadrature()`, `PetscFECreateDefault()`, `PetscDTGaussTensorQuadrature()`, `PetscDTSimplexQuadrature()`, `PetscDTTensorQuadratureCreate()`

# External Links
$(_doc_external("DM/PetscDTCreateQuadratureByCell"))
"""
function PetscDTCreateQuadratureByCell(petsclib::PetscLibType, ct::DMPolytopeType, qorder::PetscInt, qtype::PetscDTSimplexQuadratureType) end

@for_petsc function PetscDTCreateQuadratureByCell(petsclib::$UnionPetscLib, ct::DMPolytopeType, qorder::$PetscInt, qtype::PetscDTSimplexQuadratureType )
	q_ = Ref{PetscQuadrature}()
	fq_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:PetscDTCreateQuadratureByCell, $petsc_library),
               PetscErrorCode,
               (DMPolytopeType, $PetscInt, PetscDTSimplexQuadratureType, Ptr{PetscQuadrature}, Ptr{PetscQuadrature}),
               ct, qorder, qtype, q_, fq_,
              )

	q = q_[]
	fq = fq_[]

	return q,fq
end 

"""
	q::PetscQuadrature,fq::PetscQuadrature = PetscDTCreateDefaultQuadrature(petsclib::PetscLibType,ct::DMPolytopeType, qorder::PetscInt) 
Create default quadrature for a given cell

Not collective

Input Parameters:
- `ct`     - The integration domain
- `qorder` - The desired quadrature order

Output Parameters:
- `q`  - The cell quadrature
- `fq` - The face quadrature

Level: developer

-seealso: `PetscDTCreateQuadratureByCell()`, `PetscFECreateDefault()`, `PetscDTGaussTensorQuadrature()`, `PetscDTSimplexQuadrature()`, `PetscDTTensorQuadratureCreate()`

# External Links
$(_doc_external("DM/PetscDTCreateDefaultQuadrature"))
"""
function PetscDTCreateDefaultQuadrature(petsclib::PetscLibType, ct::DMPolytopeType, qorder::PetscInt) end

@for_petsc function PetscDTCreateDefaultQuadrature(petsclib::$UnionPetscLib, ct::DMPolytopeType, qorder::$PetscInt )
	q_ = Ref{PetscQuadrature}()
	fq_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:PetscDTCreateDefaultQuadrature, $petsc_library),
               PetscErrorCode,
               (DMPolytopeType, $PetscInt, Ptr{PetscQuadrature}, Ptr{PetscQuadrature}),
               ct, qorder, q_, fq_,
              )

	q = q_[]
	fq = fq_[]

	return q,fq
end 

"""
	PetscHDF5IntCast(petsclib::PetscLibType,a::PetscInt, b::hCsize_t) 

# External Links
$(_doc_external("Viewer/PetscHDF5IntCast"))
"""
function PetscHDF5IntCast(petsclib::PetscLibType, a::PetscInt, b::hCsize_t) end

@for_petsc function PetscHDF5IntCast(petsclib::$UnionPetscLib, a::$PetscInt, b::hCsize_t )

    @chk ccall(
               (:PetscHDF5IntCast, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{hCsize_t}),
               a, b,
              )


	return nothing
end 

"""
	PetscStrtolower(petsclib::PetscLibType,a::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrtolower"))
"""
function PetscStrtolower(petsclib::PetscLibType, a::Vector{Cchar}) end

@for_petsc function PetscStrtolower(petsclib::$UnionPetscLib, a::Vector{Cchar} )

    @chk ccall(
               (:PetscStrtolower, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               a,
              )


	return nothing
end 

"""
	PetscStrtoupper(petsclib::PetscLibType,a::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrtoupper"))
"""
function PetscStrtoupper(petsclib::PetscLibType, a::Vector{Cchar}) end

@for_petsc function PetscStrtoupper(petsclib::$UnionPetscLib, a::Vector{Cchar} )

    @chk ccall(
               (:PetscStrtoupper, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar},),
               a,
              )


	return nothing
end 

"""
	PetscStrlen(petsclib::PetscLibType,s::Vector{Cchar}, len::Csize_t) 

# External Links
$(_doc_external("Sys/PetscStrlen"))
"""
function PetscStrlen(petsclib::PetscLibType, s::Vector{Cchar}, len::Csize_t) end

@for_petsc function PetscStrlen(petsclib::$UnionPetscLib, s::Vector{Cchar}, len::Csize_t )

    @chk ccall(
               (:PetscStrlen, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Csize_t}),
               s, len,
              )


	return nothing
end 

"""
	PetscStrallocpy(petsclib::PetscLibType,s::Vector{Cchar}, t::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrallocpy"))
"""
function PetscStrallocpy(petsclib::PetscLibType, s::Vector{Cchar}, t::Vector{Cchar}) end

@for_petsc function PetscStrallocpy(petsclib::$UnionPetscLib, s::Vector{Cchar}, t::Vector{Cchar} )
	t_ = Ref(pointer(t))

    @chk ccall(
               (:PetscStrallocpy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Ptr{Cchar}}),
               s, t_,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscStrcmp(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrcmp"))
"""
function PetscStrcmp(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}) end

@for_petsc function PetscStrcmp(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrcmp, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               a, b, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscStrncpy(petsclib::PetscLibType,s::Vector{Cchar}, t::Vector{Cchar}, n::Csize_t) 

# External Links
$(_doc_external("Sys/PetscStrncpy"))
"""
function PetscStrncpy(petsclib::PetscLibType, s::Vector{Cchar}, t::Vector{Cchar}, n::Csize_t) end

@for_petsc function PetscStrncpy(petsclib::$UnionPetscLib, s::Vector{Cchar}, t::Vector{Cchar}, n::Csize_t )

    @chk ccall(
               (:PetscStrncpy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               s, t, n,
              )


	return nothing
end 

"""
	PetscStrlcat(petsclib::PetscLibType,s::Vector{Cchar}, t::Vector{Cchar}, n::Csize_t) 

# External Links
$(_doc_external("Sys/PetscStrlcat"))
"""
function PetscStrlcat(petsclib::PetscLibType, s::Vector{Cchar}, t::Vector{Cchar}, n::Csize_t) end

@for_petsc function PetscStrlcat(petsclib::$UnionPetscLib, s::Vector{Cchar}, t::Vector{Cchar}, n::Csize_t )

    @chk ccall(
               (:PetscStrlcat, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Csize_t),
               s, t, n,
              )


	return nothing
end 

"""
	t::PetscBool = PetscStrncmp(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}, n::Csize_t) 

# External Links
$(_doc_external("Sys/PetscStrncmp"))
"""
function PetscStrncmp(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}, n::Csize_t) end

@for_petsc function PetscStrncmp(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar}, n::Csize_t )
	t_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrncmp, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               a, b, n, t_,
              )

	t = t_[]

	return t
end 

"""
	PetscStrrstr(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}, tmp::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrrstr"))
"""
function PetscStrrstr(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}, tmp::Vector{Cchar}) end

@for_petsc function PetscStrrstr(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar}, tmp::Vector{Cchar} )
	tmp_ = Ref(pointer(tmp))

    @chk ccall(
               (:PetscStrrstr, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
               a, b, tmp_,
              )


	return nothing
end 

"""
	PetscStrstr(petsclib::PetscLibType,haystack::Vector{Cchar}, needle::Vector{Cchar}, tmp::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrstr"))
"""
function PetscStrstr(petsclib::PetscLibType, haystack::Vector{Cchar}, needle::Vector{Cchar}, tmp::Vector{Cchar}) end

@for_petsc function PetscStrstr(petsclib::$UnionPetscLib, haystack::Vector{Cchar}, needle::Vector{Cchar}, tmp::Vector{Cchar} )
	tmp_ = Ref(pointer(tmp))

    @chk ccall(
               (:PetscStrstr, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
               haystack, needle, tmp_,
              )


	return nothing
end 

"""
	t::PetscBool = PetscStrgrt(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrgrt"))
"""
function PetscStrgrt(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}) end

@for_petsc function PetscStrgrt(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar} )
	t_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrgrt, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               a, b, t_,
              )

	t = t_[]

	return t
end 

"""
	PetscStrchr(petsclib::PetscLibType,a::Vector{Cchar}, b::Cchar, c::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrchr"))
"""
function PetscStrchr(petsclib::PetscLibType, a::Vector{Cchar}, b::Cchar, c::Vector{Cchar}) end

@for_petsc function PetscStrchr(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Cchar, c::Vector{Cchar} )
	c_ = Ref(pointer(c))

    @chk ccall(
               (:PetscStrchr, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{Ptr{Cchar}}),
               a, b, c_,
              )


	return nothing
end 

"""
	PetscStrrchr(petsclib::PetscLibType,a::Vector{Cchar}, b::Cchar, c::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrrchr"))
"""
function PetscStrrchr(petsclib::PetscLibType, a::Vector{Cchar}, b::Cchar, c::Vector{Cchar}) end

@for_petsc function PetscStrrchr(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Cchar, c::Vector{Cchar} )
	c_ = Ref(pointer(c))

    @chk ccall(
               (:PetscStrrchr, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{Ptr{Cchar}}),
               a, b, c_,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscStrendswith(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrendswith"))
"""
function PetscStrendswith(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}) end

@for_petsc function PetscStrendswith(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrendswith, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               a, b, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PetscStrbeginswith(petsclib::PetscLibType,a::Vector{Cchar}, b::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscStrbeginswith"))
"""
function PetscStrbeginswith(petsclib::PetscLibType, a::Vector{Cchar}, b::Vector{Cchar}) end

@for_petsc function PetscStrbeginswith(petsclib::$UnionPetscLib, a::Vector{Cchar}, b::Vector{Cchar} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscStrbeginswith, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               a, b, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscMemmove(petsclib::PetscLibType,a::Cvoid, b::Cvoid, n::Csize_t) 

# External Links
$(_doc_external("Sys/PetscMemmove"))
"""
function PetscMemmove(petsclib::PetscLibType, a::Cvoid, b::Cvoid, n::Csize_t) end

@for_petsc function PetscMemmove(petsclib::$UnionPetscLib, a::Cvoid, b::Cvoid, n::Csize_t )

    @chk ccall(
               (:PetscMemmove, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
               a, b, n,
              )


	return nothing
end 

"""
	PetscMemcpy(petsclib::PetscLibType,a::Cvoid, b::Cvoid, n::Csize_t) 

# External Links
$(_doc_external("Sys/PetscMemcpy"))
"""
function PetscMemcpy(petsclib::PetscLibType, a::Cvoid, b::Cvoid, n::Csize_t) end

@for_petsc function PetscMemcpy(petsclib::$UnionPetscLib, a::Cvoid, b::Cvoid, n::Csize_t )

    @chk ccall(
               (:PetscMemcpy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
               a, b, n,
              )


	return nothing
end 

"""
	PetscMemzero(petsclib::PetscLibType,a::Cvoid, n::Csize_t) 

# External Links
$(_doc_external("Sys/PetscMemzero"))
"""
function PetscMemzero(petsclib::PetscLibType, a::Cvoid, n::Csize_t) end

@for_petsc function PetscMemzero(petsclib::$UnionPetscLib, a::Cvoid, n::Csize_t )

    @chk ccall(
               (:PetscMemzero, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, Csize_t),
               a, n,
              )


	return nothing
end 

"""
	PetscBTMemzero(petsclib::PetscLibType,m::PetscCount, array::PetscBT) 

# External Links
$(_doc_external("Sys/PetscBTMemzero"))
"""
function PetscBTMemzero(petsclib::PetscLibType, m::PetscCount, array::PetscBT) end

@for_petsc function PetscBTMemzero(petsclib::$UnionPetscLib, m::PetscCount, array::PetscBT )

    @chk ccall(
               (:PetscBTMemzero, $petsc_library),
               PetscErrorCode,
               (PetscCount, PetscBT),
               m, array,
              )


	return nothing
end 

"""
	PetscBTDestroy(petsclib::PetscLibType,array::PetscBT) 

# External Links
$(_doc_external("Sys/PetscBTDestroy"))
"""
function PetscBTDestroy(petsclib::PetscLibType, array::PetscBT) end

@for_petsc function PetscBTDestroy(petsclib::$UnionPetscLib, array::PetscBT )

    @chk ccall(
               (:PetscBTDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBT},),
               array,
              )


	return nothing
end 

"""
	array::PetscBT = PetscBTCreate(petsclib::PetscLibType,m::PetscCount) 

# External Links
$(_doc_external("Sys/PetscBTCreate"))
"""
function PetscBTCreate(petsclib::PetscLibType, m::PetscCount) end

@for_petsc function PetscBTCreate(petsclib::$UnionPetscLib, m::PetscCount )
	array_ = Ref{PetscBT}()

    @chk ccall(
               (:PetscBTCreate, $petsc_library),
               PetscErrorCode,
               (PetscCount, Ptr{PetscBT}),
               m, array_,
              )

	array = array_[]

	return array
end 

"""
	PetscBTCopy(petsclib::PetscLibType,dest::PetscBT, m::PetscCount, source::PetscBT) 

# External Links
$(_doc_external("Sys/PetscBTCopy"))
"""
function PetscBTCopy(petsclib::PetscLibType, dest::PetscBT, m::PetscCount, source::PetscBT) end

@for_petsc function PetscBTCopy(petsclib::$UnionPetscLib, dest::PetscBT, m::PetscCount, source::PetscBT )

    @chk ccall(
               (:PetscBTCopy, $petsc_library),
               PetscErrorCode,
               (PetscBT, PetscCount, PetscBT),
               dest, m, source,
              )


	return nothing
end 

"""
	PetscBTSet(petsclib::PetscLibType,array::PetscBT, index::PetscCount) 

# External Links
$(_doc_external("Sys/PetscBTSet"))
"""
function PetscBTSet(petsclib::PetscLibType, array::PetscBT, index::PetscCount) end

@for_petsc function PetscBTSet(petsclib::$UnionPetscLib, array::PetscBT, index::PetscCount )

    @chk ccall(
               (:PetscBTSet, $petsc_library),
               PetscErrorCode,
               (PetscBT, PetscCount),
               array, index,
              )


	return nothing
end 

"""
	PetscBTNegate(petsclib::PetscLibType,array::PetscBT, index::PetscCount) 

# External Links
$(_doc_external("Sys/PetscBTNegate"))
"""
function PetscBTNegate(petsclib::PetscLibType, array::PetscBT, index::PetscCount) end

@for_petsc function PetscBTNegate(petsclib::$UnionPetscLib, array::PetscBT, index::PetscCount )

    @chk ccall(
               (:PetscBTNegate, $petsc_library),
               PetscErrorCode,
               (PetscBT, PetscCount),
               array, index,
              )


	return nothing
end 

"""
	PetscBTClear(petsclib::PetscLibType,array::PetscBT, index::PetscCount) 

# External Links
$(_doc_external("Sys/PetscBTClear"))
"""
function PetscBTClear(petsclib::PetscLibType, array::PetscBT, index::PetscCount) end

@for_petsc function PetscBTClear(petsclib::$UnionPetscLib, array::PetscBT, index::PetscCount )

    @chk ccall(
               (:PetscBTClear, $petsc_library),
               PetscErrorCode,
               (PetscBT, PetscCount),
               array, index,
              )


	return nothing
end 

"""
	PetscLogGpuFlops(petsclib::PetscLibType,n::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogGpuFlops"))
"""
function PetscLogGpuFlops(petsclib::PetscLibType, n::PetscLogDouble) end

@for_petsc function PetscLogGpuFlops(petsclib::$UnionPetscLib, n::PetscLogDouble )

    @chk ccall(
               (:PetscLogGpuFlops, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               n,
              )


	return nothing
end 

"""
	PetscLogGpuTimeAdd(petsclib::PetscLibType,t::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogGpuTimeAdd"))
"""
function PetscLogGpuTimeAdd(petsclib::PetscLibType, t::PetscLogDouble) end

@for_petsc function PetscLogGpuTimeAdd(petsclib::$UnionPetscLib, t::PetscLogDouble )

    @chk ccall(
               (:PetscLogGpuTimeAdd, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               t,
              )


	return nothing
end 

"""
	PetscLogCpuToGpu(petsclib::PetscLibType,size::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogCpuToGpu"))
"""
function PetscLogCpuToGpu(petsclib::PetscLibType, size::PetscLogDouble) end

@for_petsc function PetscLogCpuToGpu(petsclib::$UnionPetscLib, size::PetscLogDouble )

    @chk ccall(
               (:PetscLogCpuToGpu, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               size,
              )


	return nothing
end 

"""
	PetscLogGpuToCpu(petsclib::PetscLibType,size::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogGpuToCpu"))
"""
function PetscLogGpuToCpu(petsclib::PetscLibType, size::PetscLogDouble) end

@for_petsc function PetscLogGpuToCpu(petsclib::$UnionPetscLib, size::PetscLogDouble )

    @chk ccall(
               (:PetscLogGpuToCpu, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               size,
              )


	return nothing
end 

"""
	PetscLogCpuToGpuScalar(petsclib::PetscLibType,size::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogCpuToGpuScalar"))
"""
function PetscLogCpuToGpuScalar(petsclib::PetscLibType, size::PetscLogDouble) end

@for_petsc function PetscLogCpuToGpuScalar(petsclib::$UnionPetscLib, size::PetscLogDouble )

    @chk ccall(
               (:PetscLogCpuToGpuScalar, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               size,
              )


	return nothing
end 

"""
	PetscLogGpuToCpuScalar(petsclib::PetscLibType,size::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscLogGpuToCpuScalar"))
"""
function PetscLogGpuToCpuScalar(petsclib::PetscLibType, size::PetscLogDouble) end

@for_petsc function PetscLogGpuToCpuScalar(petsclib::$UnionPetscLib, size::PetscLogDouble )

    @chk ccall(
               (:PetscLogGpuToCpuScalar, $petsc_library),
               PetscErrorCode,
               (PetscLogDouble,),
               size,
              )


	return nothing
end 

"""
	b::PetscInt = PetscIntCast(petsclib::PetscLibType,a::MPIU_Count) 

# External Links
$(_doc_external("Sys/PetscIntCast"))
"""
function PetscIntCast(petsclib::PetscLibType, a::MPIU_Count) end

@for_petsc function PetscIntCast(petsclib::$UnionPetscLib, a::MPIU_Count )
	b_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscIntCast, $petsc_library),
               PetscErrorCode,
               (MPIU_Count, Ptr{$PetscInt}),
               a, b_,
              )

	b = b_[]

	return b
end 

"""
	PetscBLASIntCast(petsclib::PetscLibType,a::MPIU_Count, b::PetscBLASInt) 

# External Links
$(_doc_external("Sys/PetscBLASIntCast"))
"""
function PetscBLASIntCast(petsclib::PetscLibType, a::MPIU_Count, b::PetscBLASInt) end

@for_petsc function PetscBLASIntCast(petsclib::$UnionPetscLib, a::MPIU_Count, b::PetscBLASInt )

    @chk ccall(
               (:PetscBLASIntCast, $petsc_library),
               PetscErrorCode,
               (MPIU_Count, Ptr{PetscBLASInt}),
               a, b,
              )


	return nothing
end 

"""
	PetscCuBLASIntCast(petsclib::PetscLibType,a::MPIU_Count, b::PetscCuBLASInt) 

# External Links
$(_doc_external("Sys/PetscCuBLASIntCast"))
"""
function PetscCuBLASIntCast(petsclib::PetscLibType, a::MPIU_Count, b::PetscCuBLASInt) end

@for_petsc function PetscCuBLASIntCast(petsclib::$UnionPetscLib, a::MPIU_Count, b::PetscCuBLASInt )

    @chk ccall(
               (:PetscCuBLASIntCast, $petsc_library),
               PetscErrorCode,
               (MPIU_Count, Ptr{PetscCuBLASInt}),
               a, b,
              )


	return nothing
end 

"""
	PetscHipBLASIntCast(petsclib::PetscLibType,a::MPIU_Count, b::PetscHipBLASInt) 

# External Links
$(_doc_external("Sys/PetscHipBLASIntCast"))
"""
function PetscHipBLASIntCast(petsclib::PetscLibType, a::MPIU_Count, b::PetscHipBLASInt) end

@for_petsc function PetscHipBLASIntCast(petsclib::$UnionPetscLib, a::MPIU_Count, b::PetscHipBLASInt )

    @chk ccall(
               (:PetscHipBLASIntCast, $petsc_library),
               PetscErrorCode,
               (MPIU_Count, Ptr{PetscHipBLASInt}),
               a, b,
              )


	return nothing
end 

"""
	PetscMPIIntCast(petsclib::PetscLibType,a::MPIU_Count, b::PetscMPIInt) 

# External Links
$(_doc_external("Sys/PetscMPIIntCast"))
"""
function PetscMPIIntCast(petsclib::PetscLibType, a::MPIU_Count, b::PetscMPIInt) end

@for_petsc function PetscMPIIntCast(petsclib::$UnionPetscLib, a::MPIU_Count, b::PetscMPIInt )

    @chk ccall(
               (:PetscMPIIntCast, $petsc_library),
               PetscErrorCode,
               (MPIU_Count, Ptr{PetscMPIInt}),
               a, b,
              )


	return nothing
end 

"""
	PetscCIntCast(petsclib::PetscLibType,a::MPIU_Count, b::Cint) 

# External Links
$(_doc_external("Sys/PetscCIntCast"))
"""
function PetscCIntCast(petsclib::PetscLibType, a::MPIU_Count, b::Cint) end

@for_petsc function PetscCIntCast(petsclib::$UnionPetscLib, a::MPIU_Count, b::Cint )

    @chk ccall(
               (:PetscCIntCast, $petsc_library),
               PetscErrorCode,
               (MPIU_Count, Ptr{Cint}),
               a, b,
              )


	return nothing
end 

"""
	result::PetscInt = PetscIntMultError(petsclib::PetscLibType,a::PetscInt, b::PetscInt) 

# External Links
$(_doc_external("Sys/PetscIntMultError"))
"""
function PetscIntMultError(petsclib::PetscLibType, a::PetscInt, b::PetscInt) end

@for_petsc function PetscIntMultError(petsclib::$UnionPetscLib, a::$PetscInt, b::$PetscInt )
	result_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscIntMultError, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}),
               a, b, result_,
              )

	result = result_[]

	return result
end 

"""
	result::PetscInt = PetscIntSumError(petsclib::PetscLibType,a::PetscInt, b::PetscInt) 

# External Links
$(_doc_external("Sys/PetscIntSumError"))
"""
function PetscIntSumError(petsclib::PetscLibType, a::PetscInt, b::PetscInt) end

@for_petsc function PetscIntSumError(petsclib::$UnionPetscLib, a::$PetscInt, b::$PetscInt )
	result_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscIntSumError, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}),
               a, b, result_,
              )

	result = result_[]

	return result
end 

"""
	set::PetscBool = PetscCitationsRegister(petsclib::PetscLibType,cit::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscCitationsRegister"))
"""
function PetscCitationsRegister(petsclib::PetscLibType, cit::Vector{Cchar}) end

@for_petsc function PetscCitationsRegister(petsclib::$UnionPetscLib, cit::Vector{Cchar} )
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscCitationsRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscBool}),
               cit, set_,
              )

	set = set_[]

	return set
end 

"""
	factorial::PetscReal = PetscDTFactorial(petsclib::PetscLibType,n::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTFactorial"))
"""
function PetscDTFactorial(petsclib::PetscLibType, n::PetscInt) end

@for_petsc function PetscDTFactorial(petsclib::$UnionPetscLib, n::$PetscInt )
	factorial_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTFactorial, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscReal}),
               n, factorial_,
              )

	factorial = factorial_[]

	return factorial
end 

"""
	factorial::PetscInt = PetscDTFactorialInt(petsclib::PetscLibType,n::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTFactorialInt"))
"""
function PetscDTFactorialInt(petsclib::PetscLibType, n::PetscInt) end

@for_petsc function PetscDTFactorialInt(petsclib::$UnionPetscLib, n::$PetscInt )
	factorial_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTFactorialInt, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}),
               n, factorial_,
              )

	factorial = factorial_[]

	return factorial
end 

"""
	binomial::PetscReal = PetscDTBinomial(petsclib::PetscLibType,n::PetscInt, k::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTBinomial"))
"""
function PetscDTBinomial(petsclib::PetscLibType, n::PetscInt, k::PetscInt) end

@for_petsc function PetscDTBinomial(petsclib::$UnionPetscLib, n::$PetscInt, k::$PetscInt )
	binomial_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDTBinomial, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscReal}),
               n, k, binomial_,
              )

	binomial = binomial_[]

	return binomial
end 

"""
	binomial::PetscInt = PetscDTBinomialInt(petsclib::PetscLibType,n::PetscInt, k::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTBinomialInt"))
"""
function PetscDTBinomialInt(petsclib::PetscLibType, n::PetscInt, k::PetscInt) end

@for_petsc function PetscDTBinomialInt(petsclib::$UnionPetscLib, n::$PetscInt, k::$PetscInt )
	binomial_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTBinomialInt, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}),
               n, k, binomial_,
              )

	binomial = binomial_[]

	return binomial
end 

"""
	perm::PetscInt,isOdd::PetscBool = PetscDTEnumPerm(petsclib::PetscLibType,n::PetscInt, k::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTEnumPerm"))
"""
function PetscDTEnumPerm(petsclib::PetscLibType, n::PetscInt, k::PetscInt) end

@for_petsc function PetscDTEnumPerm(petsclib::$UnionPetscLib, n::$PetscInt, k::$PetscInt )
	perm_ = Ref{$PetscInt}()
	isOdd_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDTEnumPerm, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, k, perm_, isOdd_,
              )

	perm = perm_[]
	isOdd = isOdd_[]

	return perm,isOdd
end 

"""
	perm::PetscInt,k::PetscInt,isOdd::PetscBool = PetscDTPermIndex(petsclib::PetscLibType,n::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTPermIndex"))
"""
function PetscDTPermIndex(petsclib::PetscLibType, n::PetscInt) end

@for_petsc function PetscDTPermIndex(petsclib::$UnionPetscLib, n::$PetscInt )
	perm_ = Ref{$PetscInt}()
	k_ = Ref{$PetscInt}()
	isOdd_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDTPermIndex, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, perm_, k_, isOdd_,
              )

	perm = perm_[]
	k = k_[]
	isOdd = isOdd_[]

	return perm,k,isOdd
end 

"""
	subset::PetscInt = PetscDTEnumSubset(petsclib::PetscLibType,n::PetscInt, k::PetscInt, j::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTEnumSubset"))
"""
function PetscDTEnumSubset(petsclib::PetscLibType, n::PetscInt, k::PetscInt, j::PetscInt) end

@for_petsc function PetscDTEnumSubset(petsclib::$UnionPetscLib, n::$PetscInt, k::$PetscInt, j::$PetscInt )
	subset_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTEnumSubset, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               n, k, j, subset_,
              )

	subset = subset_[]

	return subset
end 

"""
	subset::PetscInt,index::PetscInt = PetscDTSubsetIndex(petsclib::PetscLibType,n::PetscInt, k::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTSubsetIndex"))
"""
function PetscDTSubsetIndex(petsclib::PetscLibType, n::PetscInt, k::PetscInt) end

@for_petsc function PetscDTSubsetIndex(petsclib::$UnionPetscLib, n::$PetscInt, k::$PetscInt )
	subset_ = Ref{$PetscInt}()
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDTSubsetIndex, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               n, k, subset_, index_,
              )

	subset = subset_[]
	index = index_[]

	return subset,index
end 

"""
	perm::PetscInt,isOdd::PetscBool = PetscDTEnumSplit(petsclib::PetscLibType,n::PetscInt, k::PetscInt, j::PetscInt) 

# External Links
$(_doc_external("DM/PetscDTEnumSplit"))
"""
function PetscDTEnumSplit(petsclib::PetscLibType, n::PetscInt, k::PetscInt, j::PetscInt) end

@for_petsc function PetscDTEnumSplit(petsclib::$UnionPetscLib, n::$PetscInt, k::$PetscInt, j::$PetscInt )
	perm_ = Ref{$PetscInt}()
	isOdd_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDTEnumSplit, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               n, k, j, perm_, isOdd_,
              )

	perm = perm_[]
	isOdd = isOdd_[]

	return perm,isOdd
end 

"""
	PetscTime(petsclib::PetscLibType,v::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscTime"))
"""
function PetscTime(petsclib::PetscLibType, v::PetscLogDouble) end

@for_petsc function PetscTime(petsclib::$UnionPetscLib, v::PetscLogDouble )

    @chk ccall(
               (:PetscTime, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               v,
              )


	return nothing
end 

"""
	PetscTimeSubtract(petsclib::PetscLibType,v::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscTimeSubtract"))
"""
function PetscTimeSubtract(petsclib::PetscLibType, v::PetscLogDouble) end

@for_petsc function PetscTimeSubtract(petsclib::$UnionPetscLib, v::PetscLogDouble )

    @chk ccall(
               (:PetscTimeSubtract, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               v,
              )


	return nothing
end 

"""
	PetscTimeAdd(petsclib::PetscLibType,v::PetscLogDouble) 

# External Links
$(_doc_external("Sys/PetscTimeAdd"))
"""
function PetscTimeAdd(petsclib::PetscLibType, v::PetscLogDouble) end

@for_petsc function PetscTimeAdd(petsclib::$UnionPetscLib, v::PetscLogDouble )

    @chk ccall(
               (:PetscTimeAdd, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogDouble},),
               v,
              )


	return nothing
end 

"""
	PetscOptionsBegin(petsclib::PetscLibType,comm::MPI_Comm, prefix::Vector{Cchar}, mess::Vector{Cchar}, sec::Vector{Cchar}) 

# External Links
$(_doc_external("Sys/PetscOptionsBegin"))
"""
function PetscOptionsBegin(petsclib::PetscLibType, comm::MPI_Comm, prefix::Vector{Cchar}, mess::Vector{Cchar}, sec::Vector{Cchar}) end

@for_petsc function PetscOptionsBegin(petsclib::$UnionPetscLib, comm::MPI_Comm, prefix::Vector{Cchar}, mess::Vector{Cchar}, sec::Vector{Cchar} )

    @chk ccall(
               (:PetscOptionsBegin, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
               comm, prefix, mess, sec,
              )


	return nothing
end 

"""
	PetscOptionsEnd(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscOptionsEnd"))
"""
function PetscOptionsEnd(petsclib::PetscLibType) end

@for_petsc function PetscOptionsEnd(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscOptionsEnd, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	value::PetscBool,set::PetscBool = PetscOptionsBool(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscBool) 

# External Links
$(_doc_external("Sys/PetscOptionsBool"))
"""
function PetscOptionsBool(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscBool) end

@for_petsc function PetscOptionsBool(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscBool )
	value_ = Ref{PetscBool}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsBool, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, PetscBool, Ptr{PetscBool}, Ptr{PetscBool}),
               opt, text, man, currentvalue, value_, set_,
              )

	value = value_[]
	set = set_[]

	return value,set
end 

"""
	PetscOptionsBool3(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscBool3, value::PetscBool3, set::PetscBool3) 

# External Links
$(_doc_external("Sys/PetscOptionsBool3"))
"""
function PetscOptionsBool3(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscBool3, value::PetscBool3, set::PetscBool3) end

@for_petsc function PetscOptionsBool3(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscBool3, value::PetscBool3, set::PetscBool3 )

    @chk ccall(
               (:PetscOptionsBool3, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, PetscBool3, Ptr{PetscBool3}, Ptr{PetscBool3}),
               opt, text, man, currentvalue, value, set,
              )


	return nothing
end 

"""
	value::PetscInt,set::PetscBool = PetscOptionsInt(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscInt) 

# External Links
$(_doc_external("Sys/PetscOptionsInt"))
"""
function PetscOptionsInt(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscInt) end

@for_petsc function PetscOptionsInt(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::$PetscInt )
	value_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsInt, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               opt, text, man, currentvalue, value_, set_,
              )

	value = value_[]
	set = set_[]

	return value,set
end 

"""
	value::PetscReal,set::PetscBool = PetscOptionsReal(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscReal) 

# External Links
$(_doc_external("Sys/PetscOptionsReal"))
"""
function PetscOptionsReal(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscReal) end

@for_petsc function PetscOptionsReal(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::$PetscReal )
	value_ = Ref{$PetscReal}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsReal, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, $PetscReal, Ptr{$PetscReal}, Ptr{PetscBool}),
               opt, text, man, currentvalue, value_, set_,
              )

	value = value_[]
	set = set_[]

	return value,set
end 

"""
	value::PetscScalar,set::PetscBool = PetscOptionsScalar(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscScalar) 

# External Links
$(_doc_external("Sys/PetscOptionsScalar"))
"""
function PetscOptionsScalar(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::PetscScalar) end

@for_petsc function PetscOptionsScalar(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, currentvalue::$PetscScalar )
	value_ = Ref{$PetscScalar}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsScalar, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, $PetscScalar, Ptr{$PetscScalar}, Ptr{PetscBool}),
               opt, text, man, currentvalue, value_, set_,
              )

	value = value_[]
	set = set_[]

	return value,set
end 

"""
	n::PetscInt,set::PetscBool = PetscOptionsScalarArray(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscScalar}) 

# External Links
$(_doc_external("Sys/PetscOptionsScalarArray"))
"""
function PetscOptionsScalarArray(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscScalar}) end

@for_petsc function PetscOptionsScalarArray(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{$PetscScalar} )
	n_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsScalarArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscScalar}, Ptr{$PetscInt}, Ptr{PetscBool}),
               opt, text, man, value, n_, set_,
              )

	n = n_[]
	set = set_[]

	return n,set
end 

"""
	n::PetscInt,set::PetscBool = PetscOptionsIntArray(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscInt}) 

# External Links
$(_doc_external("Sys/PetscOptionsIntArray"))
"""
function PetscOptionsIntArray(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscInt}) end

@for_petsc function PetscOptionsIntArray(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsIntArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{PetscBool}),
               opt, text, man, value, n_, set_,
              )

	n = n_[]
	set = set_[]

	return n,set
end 

"""
	n::PetscInt,set::PetscBool = PetscOptionsRealArray(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscReal}) 

# External Links
$(_doc_external("Sys/PetscOptionsRealArray"))
"""
function PetscOptionsRealArray(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscReal}) end

@for_petsc function PetscOptionsRealArray(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{$PetscReal} )
	n_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsRealArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{PetscBool}),
               opt, text, man, value, n_, set_,
              )

	n = n_[]
	set = set_[]

	return n,set
end 

"""
	n::PetscInt,set::PetscBool = PetscOptionsBoolArray(petsclib::PetscLibType,opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscBool}) 

# External Links
$(_doc_external("Sys/PetscOptionsBoolArray"))
"""
function PetscOptionsBoolArray(petsclib::PetscLibType, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscBool}) end

@for_petsc function PetscOptionsBoolArray(petsclib::$UnionPetscLib, opt::Vector{Cchar}, text::Vector{Cchar}, man::Vector{Cchar}, value::Vector{PetscBool} )
	n_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsBoolArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}, Ptr{$PetscInt}, Ptr{PetscBool}),
               opt, text, man, value, n_, set_,
              )

	n = n_[]
	set = set_[]

	return n,set
end 

