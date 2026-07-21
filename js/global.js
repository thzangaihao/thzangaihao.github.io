/* ============================================================
   文件名: global.js
   描述: 全站通用脚本 (终极修复版 - 解决语法崩溃与重复加载)
   ============================================================ */

/* ------------------------------------------------------------
   1. 全局配置 (Configuration)
   ------------------------------------------------------------ */
window.MathJax = {
    tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        displayMath: [['$$', '$$'], ['\\[', '\\]']],
        tags: 'ams'
    },
    options: { ignoreHtmlClass: 'tex2jax_ignore', processHtmlClass: 'tex2jax_process' }
};

window.__siteGlobalScriptUrl = window.__siteGlobalScriptUrl || document.currentScript.src;
function getSiteImageUrl(fileName) {
    return new URL(`../img/${fileName}`, window.__siteGlobalScriptUrl).href;
}

// 在首次绘制前恢复外观，避免页面先显示浅色背景再跳变。
(() => {
    let backgroundMode = 'image';
    let themeMode = 'system';
    try {
        backgroundMode = localStorage.getItem('site-background-mode') || 'image';
        themeMode = localStorage.getItem('site-theme-mode') || 'system';
    } catch (error) {}
    const useDark = themeMode === 'dark' || (themeMode === 'system' && window.matchMedia('(prefers-color-scheme: dark)').matches);
    document.documentElement.classList.toggle('site-dark', useDark);
    document.documentElement.classList.toggle('site-background-image', backgroundMode === 'image');
    document.documentElement.style.colorScheme = useDark ? 'dark' : 'light';
    const initialPoster = useDark ? '2.png' : '1.png';
    [initialPoster].forEach(fileName => {
        const preload = document.createElement('link');
        preload.rel = 'preload';
        preload.as = 'image';
        preload.href = getSiteImageUrl(fileName);
        document.head.appendChild(preload);
    });
})();

/* ------------------------------------------------------------
   2. 全局工具函数 (Utilities)
   ------------------------------------------------------------ */
function getRelativePath(dbPath) { return '../../' + dbPath; }

function isTruthyMeta(value) {
    return ['true', 'yes', '1', 'featured', '精选'].includes((value || '').trim().toLowerCase());
}

function featuredStarHtml(article) {
    if (!article?.featured) return '';
    return '<i class="fas fa-star featured-star" title="精选文章" aria-label="精选文章"></i>';
}

let copyFeedbackTimer = null;
let currentActiveElement = null;
let originalContentCache = "";

function resetCopyFeedback() {
    if (currentActiveElement) {
        currentActiveElement.innerHTML = originalContentCache;
        currentActiveElement = null;
        originalContentCache = "";
    }
    if (copyFeedbackTimer) { clearTimeout(copyFeedbackTimer); copyFeedbackTimer = null; }
}

function copyToClipboard(text, element, event) {
    if (event) event.preventDefault();
    if (currentActiveElement && currentActiveElement !== element) resetCopyFeedback();
    if (copyFeedbackTimer) clearTimeout(copyFeedbackTimer);

    navigator.clipboard.writeText(text).then(() => {
        if (currentActiveElement !== element) {
            currentActiveElement = element;
            originalContentCache = element.innerHTML;
        }
        element.innerHTML = '<i class="fas fa-check" style="color: #27ae60;"></i> 已复制';
        copyFeedbackTimer = setTimeout(() => resetCopyFeedback(), 1500);
    }).catch(err => console.error('复制失败:', err));
}

function copyCode(text, btnElement) {
    navigator.clipboard.writeText(text).then(() => {
        const originalHtml = btnElement.innerHTML;
        btnElement.innerHTML = '<i class="fas fa-check"></i>';
        btnElement.classList.add('copied');
        setTimeout(() => {
            btnElement.innerHTML = originalHtml;
            btnElement.classList.remove('copied');
        }, 1500);
    }).catch(err => { console.error('复制失败:', err); btnElement.innerText = 'Error'; });
}

/* ------------------------------------------------------------
   3. 核心 UI 逻辑
   ------------------------------------------------------------ */
function loadFooter() {
    const path = window.footerPath || '/footer.html';
    return fetch(path)
        .then(res => { if (!res.ok) throw new Error(`Footer not found`); return res.text(); })
        .then(data => { const el = document.getElementById('footer-placeholder'); if (el) el.innerHTML = data; })
        .catch(err => console.error('Footer Error:', err));
}

function makeTablesResponsive() {
    document.querySelectorAll('.content table').forEach(table => {
        if (table.closest('.table-outer-container')) return;

        const outerContainer = document.createElement('div');
        outerContainer.className = 'table-outer-container';
        const scrollWrapper = document.createElement('div');
        scrollWrapper.className = 'table-wrapper';
        const toggleBtn = document.createElement('button');
        toggleBtn.className = 'table-toggle-btn';
        toggleBtn.innerHTML = '<i class="fas fa-align-left"></i> 自动换行';
        toggleBtn.title = '点击切换表格换行模式';
        
        toggleBtn.onclick = () => {
            table.classList.toggle('force-wrap');
            const isWrapped = table.classList.contains('force-wrap');
            if (isWrapped) {
                toggleBtn.innerHTML = '<i class="fas fa-compress"></i> 恢复不换行';
                toggleBtn.classList.add('active');
            } else {
                toggleBtn.innerHTML = '<i class="fas fa-align-left"></i> 自动换行';
                toggleBtn.classList.remove('active');
            }
        };

        const caption = table.querySelector('caption');
        let externalCaption = null;
        if (caption) {
            externalCaption = document.createElement('div');
            externalCaption.className = 'table-bottom-caption';
            externalCaption.innerHTML = caption.innerHTML;
            table.removeChild(caption); 
        }

        table.parentNode.insertBefore(outerContainer, table);
        outerContainer.appendChild(toggleBtn);
        outerContainer.appendChild(scrollWrapper);
        scrollWrapper.appendChild(table);
        if (externalCaption) {
            outerContainer.parentNode.insertBefore(externalCaption, outerContainer.nextSibling);
        }
    });
}

function autoFillArticleInfo() {
    const container = document.querySelector('.article-container');
    if (!container) return;
    const getMeta = (n) => document.querySelector(`meta[name="${n}"]`)?.content;
    const [collId, collTitle, pPath, pTitle, date] = [getMeta('collection'), getMeta('collection-title'), getMeta('parent-path'), getMeta('parent-title'), getMeta('date')];
    const isFeatured = isTruthyMeta(getMeta('featured'));

    const h1 = container.querySelector('h1');
    if (h1) {
        if (collId && collTitle && pPath && pTitle && date) {
            h1.insertAdjacentHTML('beforebegin', `<a href="${pPath}#${collId}" class="back-link"><i class="fas fa-arrow-left"></i> 返回 ${pTitle} / ${collTitle}</a>`);
            h1.insertAdjacentHTML('afterend', `<div class="meta"><span><i class="far fa-calendar"></i> ${date}</span><span style="margin-left:20px"><i class="fas fa-folder"></i> ${collTitle}</span></div>`);
        }

        if (isFeatured && !container.querySelector('.featured-notice')) {
            h1.insertAdjacentHTML('afterend', `
                <div class="featured-notice">
                    <i class="fas fa-star"></i>
                    <span>精选文章</span>
                </div>
            `);
        }
    }
}

/* ------------------------------------------------------------
   4. 导航逻辑
   ------------------------------------------------------------ */
function showCollection(collectionId, btnElement) {
    document.querySelectorAll('.article-list').forEach(el => el.classList.remove('active'));
    const target = document.getElementById(collectionId + '-section');
    if (target) target.classList.add('active');
    document.querySelectorAll('.collection-btn').forEach(el => el.classList.remove('active'));
    if (btnElement) btnElement.classList.add('active');
}

function handleHashNavigation() {
    const hash = window.location.hash.substring(1);
    if (!hash) return;
    const targetBtn = document.querySelector(`button[onclick*="'${hash}'"]`);
    if (targetBtn) showCollection(hash, targetBtn);
}

/* ------------------------------------------------------------
   5. 内容渲染 (卡片与最新文章)
   ------------------------------------------------------------ */
function renderArticleCards() {
    if (!window.ARTICLE_DATABASE || !Array.isArray(window.ARTICLE_DATABASE)) return;
    window.ARTICLE_DATABASE.forEach(article => {
        const container = document.getElementById(article.collection + '-section');
        if (!container) return;
        const html = `
            <a href="${getRelativePath(article.path)}" class="article-item${article.featured ? ' is-featured' : ''}" rel="noopener noreferrer">
                <div class="article-info">
                    <h4>${featuredStarHtml(article)}${article.title}</h4>
                    <p>${article.summary}</p>
                </div>
                <span class="article-date">${article.date}</span>
            </a>`;
        container.insertAdjacentHTML('beforeend', html);
    });
}

function renderLatestArticles(searchTerm = '') {
    const container = document.getElementById('latest-articles-list');
    const titleText = document.getElementById('section-title');
    if (!container) return;

    if (!window.ARTICLE_DATABASE) { 
        container.innerHTML = '<p style="text-align:center;color:#999">暂无数据</p>'; 
        return; 
    }
    
    if (titleText) {
        titleText.textContent = searchTerm ? '搜索结果' : '最新动态';
    }

    const lowerTerm = searchTerm.toLowerCase();
    
    // 过滤逻辑
    let latest = [...window.ARTICLE_DATABASE].filter(article => {
        return article.title.toLowerCase().includes(lowerTerm) || 
               article.summary.toLowerCase().includes(lowerTerm);
    }).sort((a, b) => new Date(b.date) - new Date(a.date));

    // 无搜索词时，限制显示最新的 10 篇文章
    if (!searchTerm) {
        latest = latest.slice(0, 10);
    }
    
    if (latest.length === 0) {
        container.innerHTML = `<div style="text-align:center; padding: 3rem 0; color:#999;"><i class="fas fa-search-minus" style="font-size:2rem; margin-bottom:10px; color:#ddd;"></i><br>未找到与 "${searchTerm}" 相关的文章</div>`;
        return;
    }

    const html = latest.map(article => {
        let sectionName = '动态'; let iconClass = 'fas fa-newspaper'; 
        if (article.path.includes('dry_lab')) { sectionName = '干实验'; iconClass = 'fas fa-code'; } 
        else if (article.path.includes('wet_lab')) { sectionName = '湿实验'; iconClass = 'fas fa-flask'; } 
        else if (article.path.includes('resources')) { sectionName = '资源站'; iconClass = 'fas fa-book'; } 
        else if (article.path.includes('about')) { sectionName = '关于我'; iconClass = 'fas fa-user'; }

        return `
            <a href="${article.path}" class="article-item${article.featured ? ' is-featured' : ''}" rel="noopener noreferrer">
                <div class="article-info">
                    <span class="article-tag"><i class="${iconClass}"></i> ${sectionName}</span>
                    <h4>${featuredStarHtml(article)}${article.title}</h4>
                    <p>${article.summary}</p>
                </div>
                <span class="article-date">${article.date}</span>
            </a>`;
    }).join('');
    container.innerHTML = html;
}

/* ------------------------------------------------------------
   6. 代码高亮与增强
   ------------------------------------------------------------ */
function loadPrismHighlighter() {
    const codeBlocks = document.querySelectorAll('pre code');
    if (codeBlocks.length === 0) return;
    const cdn = 'https://cdnjs.cloudflare.com/ajax/libs/prism/1.29.0';
    const link = document.createElement('link'); link.rel='stylesheet'; link.href=`${cdn}/themes/prism-tomorrow.min.css`; document.head.appendChild(link);
    const script = document.createElement('script'); script.src=`${cdn}/prism.min.js`; 
    script.onload = () => {
        const auto = document.createElement('script'); auto.src=`${cdn}/plugins/autoloader/prism-autoloader.min.js`;
        auto.onload = () => window.Prism.highlightAll();
        document.body.appendChild(auto);
    };
    document.body.appendChild(script);
}

function enhanceCodeBlocks() {
    const CODE_MAX_HEIGHT = 360; 

    document.querySelectorAll('pre code').forEach(code => {
        const pre = code.parentElement;
        if (pre.parentElement.classList.contains('code-wrapper')) return;

        let langName = 'TEXT';
        const langClass = Array.from(code.classList).find(c => c.startsWith('language-'));
        if (langClass) langName = langClass.replace('language-', '').toUpperCase();

        const wrapper = document.createElement('div');
        wrapper.className = 'code-wrapper';

        const header = document.createElement('div');
        header.className = 'code-header';
        header.innerHTML = `
            <span class="lang-label">${langName}</span>
            <button class="copy-code-btn" onclick="copyCode(this.parentElement.nextElementSibling.innerText, this)">
                <i class="far fa-copy"></i>
            </button>
        `;

        pre.parentNode.insertBefore(wrapper, pre);
        wrapper.appendChild(header);
        wrapper.appendChild(pre);

        const actualHeight = pre.offsetHeight;
        if (actualHeight > CODE_MAX_HEIGHT) {
            wrapper.classList.add('collapsed');
            pre.style.maxHeight = CODE_MAX_HEIGHT + 'px'; 
            
            const expandBtn = document.createElement('div');
            expandBtn.className = 'code-expand-btn';
            expandBtn.innerHTML = '<i class="fas fa-chevron-down"></i> 显示全部代码';
            
            expandBtn.onclick = () => {
                const isCollapsed = wrapper.classList.contains('collapsed');
                if (isCollapsed) {
                    wrapper.classList.remove('collapsed');
                    pre.style.maxHeight = 'none'; 
                    expandBtn.innerHTML = '<i class="fas fa-chevron-up"></i> 收起代码';
                } else {
                    wrapper.classList.add('collapsed');
                    pre.style.maxHeight = CODE_MAX_HEIGHT + 'px';
                    expandBtn.innerHTML = '<i class="fas fa-chevron-down"></i> 显示全部代码';
                    wrapper.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
                }
            };
            wrapper.appendChild(expandBtn);
        }
    });
}

/* ------------------------------------------------------------
   7. 目录生成 (TOC)
   ------------------------------------------------------------ */
function generateTOC() {
    const container = document.querySelector('.article-container');
    const content = document.querySelector('.content');
    if (!container || !content) return;

    const headers = content.querySelectorAll('h2, h3');
    if (headers.length < 2) {
        container.style.opacity = '1'; 
        return;
    }

    const wrapper = document.createElement('div');
    wrapper.className = 'article-layout-wrapper';
    container.parentNode.insertBefore(wrapper, container);
    
    const sidebar = document.createElement('aside');
    sidebar.className = 'article-sidebar';
    sidebar.innerHTML = `<div class="toc-card"><h3><i class="fas fa-list-ul"></i> 目录</h3><ul id="toc-list"></ul></div>`;
    
    wrapper.appendChild(container);
    wrapper.appendChild(sidebar);

    const list = sidebar.querySelector('#toc-list');
    
    const scrollSidebar = (link) => {
        const card = document.querySelector('.toc-card');
        if (!card) return;
        const linkTop = link.offsetTop; 
        const cardHeight = card.clientHeight;
        const scrollTop = card.scrollTop;
        if (linkTop < scrollTop + 20 || linkTop > scrollTop + cardHeight - 20) {
            card.scrollTo({ top: linkTop - 60, behavior: 'smooth' });
        }
    };

    headers.forEach((h, i) => {
        if (!h.id) h.id = 'sec-' + i;
        const li = document.createElement('li');
        if (h.tagName === 'H3') li.className = 'toc-sub-item';
        const link = document.createElement('a');
        link.href = `#${h.id}`;
        link.className = 'toc-link';
        link.dataset.target = h.id;
        link.textContent = h.textContent;
        link.onclick = (e) => {
            e.preventDefault();
            isClicking = true; 
            document.querySelectorAll('.toc-link').forEach(l => l.classList.remove('active'));
            link.classList.add('active');
            scrollSidebar(link);
            document.getElementById(h.id).scrollIntoView({ behavior: 'smooth' });
            if (window.clickTimer) clearTimeout(window.clickTimer);
            window.clickTimer = setTimeout(() => { isClicking = false; }, 1000);
        };
        li.appendChild(link);
        list.appendChild(li);
    });

    void wrapper.offsetWidth;
    wrapper.classList.add('loaded');

    let isClicking = false;
    const visibleHeaders = new Set();
    const observer = new IntersectionObserver(entries => {
        entries.forEach(entry => {
            if (entry.isIntersecting) visibleHeaders.add(entry.target.id);
            else visibleHeaders.delete(entry.target.id);
        });
        if (isClicking) return;
        let activeId = null;
        for (const header of headers) {
            if (visibleHeaders.has(header.id)) {
                activeId = header.id;
                break; 
            }
        }
        if (activeId) {
            document.querySelectorAll('.toc-link').forEach(l => {
                const isActive = l.dataset.target === activeId;
                l.classList.toggle('active', isActive);
                if (isActive) scrollSidebar(l);
            });
        }
    }, { rootMargin: '0px 0px -80% 0px', threshold: 0 });

    headers.forEach(h => observer.observe(h));
}

/* ------------------------------------------------------------
   8. 页面悬浮控制（外观设置与返回顶部）
   ------------------------------------------------------------ */
function initPageControls() {
    const backgroundStorageKey = 'site-background-mode';
    const themeStorageKey = 'site-theme-mode';
    const systemTheme = window.matchMedia('(prefers-color-scheme: dark)');
    const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');
    let backgroundMode = 'image';
    let themeMode = 'system';
    try {
        backgroundMode = localStorage.getItem(backgroundStorageKey) || 'image';
        themeMode = localStorage.getItem(themeStorageKey) || 'system';
    } catch (error) {}

    const controls = document.createElement('div');
    controls.className = 'page-fab-group';
    controls.setAttribute('aria-label', '页面显示控制');

    const appearanceBtn = document.createElement('button');
    appearanceBtn.id = 'appearance-btn';
    appearanceBtn.className = 'page-fab';
    appearanceBtn.innerHTML = '<i class="fas fa-circle-half-stroke"></i>';
    appearanceBtn.title = '外观设置';
    appearanceBtn.setAttribute('aria-label', '打开外观设置');
    appearanceBtn.setAttribute('aria-expanded', 'false');

    const panel = document.createElement('div');
    panel.id = 'appearance-panel';
    panel.className = 'appearance-panel';
    panel.hidden = true;
    panel.innerHTML = `
        <div class="appearance-heading"><span>外观</span><small>即时生效</small></div>
        <div class="appearance-setting">
            <span class="appearance-label"><i class="fas fa-circle-half-stroke"></i> 主题</span>
            <div class="appearance-options" data-setting="theme" role="group" aria-label="主题模式">
                <button type="button" data-value="light">浅色</button>
                <button type="button" data-value="dark">深色</button>
                <button type="button" data-value="system">自动</button>
            </div>
        </div>
        <div class="appearance-setting">
            <span class="appearance-label"><i class="fas fa-panorama"></i> 背景</span>
            <div class="appearance-options" data-setting="background" role="group" aria-label="背景模式">
                <button type="button" data-value="image">图片</button>
                <button type="button" data-value="solid">纯色</button>
            </div>
        </div>`;

    const topBtn = document.createElement('button');
    topBtn.id = 'back-to-top-btn';
    topBtn.className = 'page-fab';
    topBtn.innerHTML = '<i class="fas fa-arrow-up"></i>';
    topBtn.title = '回到顶部';
    topBtn.setAttribute('aria-label', '回到顶部');

    function updateOptions() {
        panel.querySelectorAll('[data-setting="theme"] button').forEach(button => {
            button.classList.toggle('active', button.dataset.value === themeMode);
        });
        panel.querySelectorAll('[data-setting="background"] button').forEach(button => {
            button.classList.toggle('active', button.dataset.value === backgroundMode);
        });
    }

    let transitionTimer = null;
    function applyAppearance(animate = false) {
        const useImage = backgroundMode === 'image';
        const useDark = themeMode === 'dark' || (themeMode === 'system' && systemTheme.matches);
        if (animate && useImage && !reducedMotion.matches) {
            document.documentElement.classList.remove('theme-transitioning');
            void document.documentElement.offsetWidth;
            document.documentElement.classList.add('theme-transitioning');
            if (transitionTimer) clearTimeout(transitionTimer);
            transitionTimer = setTimeout(() => document.documentElement.classList.remove('theme-transitioning'), 1250);
        }
        document.documentElement.classList.toggle('site-background-image', useImage);
        document.documentElement.classList.toggle('site-dark', useDark);
        document.documentElement.style.colorScheme = useDark ? 'dark' : 'light';
        appearanceBtn.classList.toggle('active', useDark);
        updateOptions();
    }

    applyAppearance();
    panel.addEventListener('click', event => {
        const option = event.target.closest('button[data-value]');
        if (!option) return;
        const setting = option.closest('[data-setting]').dataset.setting;
        if (setting === 'theme') themeMode = option.dataset.value;
        if (setting === 'background') backgroundMode = option.dataset.value;
        try {
            localStorage.setItem(themeStorageKey, themeMode);
            localStorage.setItem(backgroundStorageKey, backgroundMode);
        } catch (error) {}
        applyAppearance(true);
    });
    systemTheme.addEventListener?.('change', () => { if (themeMode === 'system') applyAppearance(true); });

    function closePanel() {
        panel.hidden = true;
        appearanceBtn.setAttribute('aria-expanded', 'false');
    }

    appearanceBtn.addEventListener('click', event => {
        event.stopPropagation();
        panel.hidden = !panel.hidden;
        appearanceBtn.setAttribute('aria-expanded', String(!panel.hidden));
    });
    document.addEventListener('click', event => { if (!controls.contains(event.target)) closePanel(); });
    document.addEventListener('keydown', event => { if (event.key === 'Escape') closePanel(); });

    controls.append(panel, appearanceBtn, topBtn);
    document.body.appendChild(controls);

    let scrollTimeout;
    window.addEventListener('scroll', () => {
        if (scrollTimeout) cancelAnimationFrame(scrollTimeout);
        scrollTimeout = requestAnimationFrame(() => {
            if (window.scrollY > 300) {
                topBtn.classList.add('show');
            } else {
                topBtn.classList.remove('show');
            }
        });
    });

    topBtn.addEventListener('click', () => {
        window.scrollTo({ top: 0, behavior: 'smooth' });
    });
}

/* ------------------------------------------------------------
   9. 自动标题编号 (Auto Numbering for H2/H3)
   ------------------------------------------------------------ */
function autoNumberHeadings() {
    const content = document.querySelector('.content');
    if (!content) return;

    // 抓取所有 h2 和 h3
    const headers = content.querySelectorAll('h2, h3');
    if (headers.length === 0) return;

    let h2Count = 0;
    let h3Count = 0;

    headers.forEach(header => {
        // 防止被重复执行导致多次编号
        if (header.classList.contains('numbered')) return;
        header.classList.add('numbered');

        if (header.tagName === 'H2') {
            h2Count++;
            h3Count = 0; // 遇到新的 H2，重置 H3 计数器
            // 插入编号 span (可以用来单独设置 CSS 样式)
            header.innerHTML = `<span class="heading-number">${h2Count}. </span>` + header.innerHTML;
        } 
        else if (header.tagName === 'H3') {
            h3Count++;
            // 如果文章一上来没有 H2 直接是 H3，为了美观不显示 0.1，直接跳过或者只显示 H3 计数
            const prefix = h2Count > 0 ? `${h2Count}.${h3Count}` : `${h3Count}`;
            header.innerHTML = `<span class="heading-number">${prefix} </span>` + header.innerHTML;
        }
    });
}

/* ------------------------------------------------------------
    初始化入口 (统一控制台)
   ------------------------------------------------------------ */
window.addEventListener('DOMContentLoaded', async () => {
    if ('scrollRestoration' in history) history.scrollRestoration = 'manual';
    initPageControls();

    await loadFooter();

    if (document.body.classList.contains('category-page')) {
        renderArticleCards(); 
        handleHashNavigation();
    }
    
    renderLatestArticles();   
    
    const searchInput = document.getElementById('article-search');
    if (searchInput) {
        searchInput.addEventListener('input', (e) => {
            renderLatestArticles(e.target.value);
        });
    }
    
    enhanceCodeBlocks();      
    loadPrismHighlighter(); 
    makeTablesResponsive();
    autoFillArticleInfo();
    autoNumberHeadings();
    generateTOC();

});
