function err = errors(rec, orig, type)

if strcmp(type, 'l2')
    err = sum(abs(rec(:) - orig(:)).^2)/sum(abs(orig(:)).^2);
elseif strcmp(type, 'SSIM')
    err = ssim(rec, orig);
elseif strcmp(type, 'PSNR')
    err = psnr(rec, orig);
end

end

